from __future__ import annotations

import unittest

from brainc._canonical import SAFE_INTEGER
from brainc.insdc_location import (
    Between,
    Complement,
    Endpoint,
    INSDCLocationError,
    Join,
    LocationLimits,
    NormalizedSegment,
    Order,
    Point,
    Remote,
    SequenceContext,
    Span,
    Within,
    normalize_segments,
    parse_location,
)


LINEAR = SequenceContext(100, "linear")
CIRCULAR = SequenceContext(100, "circular")


class INSDCLocationTests(unittest.TestCase):
    def parse(self, source: str, **kwargs):
        return parse_location(source, context=kwargs.pop("context", LINEAR), **kwargs)

    def test_points_spans_fuzz_and_zero_based_normalization(self) -> None:
        self.assertEqual(self.parse("7"), Point(7))
        location = self.parse("<2..>91")
        self.assertEqual(location, Span(Endpoint(2, "<"), Endpoint(91, ">")))
        self.assertEqual(
            normalize_segments(location, context=LINEAR),
            (NormalizedSegment(None, 1, 91, 1, "interval", "<", ">"),),
        )

    def test_join_order_and_complement_preserve_ast_child_order(self) -> None:
        joined = self.parse("join(2..4,9,complement(20..22))")
        self.assertEqual(
            joined,
            Join(
                (
                    Span(Endpoint(2), Endpoint(4)),
                    Point(9),
                    Complement(Span(Endpoint(20), Endpoint(22))),
                )
            ),
        )
        self.assertIsInstance(self.parse("order(2,4)"), Order)
        complemented = self.parse("complement(join(2..4,9..10))")
        self.assertEqual(
            normalize_segments(complemented, context=LINEAR),
            (
                NormalizedSegment(None, 8, 10, -1, "interval"),
                NormalizedSegment(None, 1, 4, -1, "interval"),
            ),
        )

    def test_between_bases_are_adjacent_or_the_explicit_circular_origin(self) -> None:
        self.assertEqual(self.parse("7^8"), Between(7, 8))
        self.assertEqual(
            normalize_segments(Between(7, 8), context=LINEAR),
            (NormalizedSegment(None, 7, 7, 1, "between"),),
        )
        origin = self.parse("100^1", context=CIRCULAR)
        self.assertEqual(
            normalize_segments(origin, context=CIRCULAR),
            (NormalizedSegment(None, 0, 0, 1, "between"),),
        )
        for source, context in (("7^9", LINEAR), ("100^1", LINEAR), ("99^1", CIRCULAR)):
            with self.subTest(source=source, topology=context.topology):
                with self.assertRaisesRegex(INSDCLocationError, "adjacent"):
                    self.parse(source, context=context)

    def test_remote_references_are_exact_and_never_implicitly_resolved(self) -> None:
        remote = self.parse("NC_000001.11:3..8")
        self.assertEqual(
            remote,
            Remote("NC_000001.11", Span(Endpoint(3), Endpoint(8))),
        )
        self.assertEqual(
            normalize_segments(remote, context=LINEAR),
            (
                NormalizedSegment(
                    "NC_000001.11",
                    2,
                    8,
                    1,
                    "interval",
                    bounds_status="unresolved",
                ),
            ),
        )
        references = {"NC_000001.11": SequenceContext(10, "linear")}
        self.assertEqual(
            normalize_segments(remote, context=LINEAR, references=references),
            (NormalizedSegment("NC_000001.11", 2, 8, 1, "interval"),),
        )
        with self.assertRaisesRegex(INSDCLocationError, "exceeds its explicit"):
            self.parse("NC_000001.11:3..11", references=references)

        bad = [
            "NC_000001:1",
            "nc_000001.11:1",
            "NC_000001.0:1",
            "NC_000001.01:1",
            "NC_000001.11:join(1,2)",
        ]
        for source in bad:
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source)

    def test_remote_circular_origin_requires_explicit_remote_topology(self) -> None:
        source = "AB000001.1:100^1"
        with self.assertRaisesRegex(INSDCLocationError, "adjacent"):
            self.parse(source)
        unresolved_adjacent = self.parse("AB000001.1:100^101")
        self.assertEqual(
            normalize_segments(unresolved_adjacent, context=LINEAR),
            (
                NormalizedSegment(
                    "AB000001.1",
                    100,
                    100,
                    1,
                    "between",
                    bounds_status="unresolved",
                ),
            ),
        )
        references = {"AB000001.1": SequenceContext(100, "circular")}
        parsed = self.parse(source, references=references)
        self.assertEqual(
            normalize_segments(parsed, context=LINEAR, references=references),
            (NormalizedSegment("AB000001.1", 0, 0, 1, "between"),),
        )

    def test_coordinates_are_ascii_positive_safe_and_bounded(self) -> None:
        invalid = [
            "0",
            "01",
            "1..02",
            str(SAFE_INTEGER + 1),
            "101",
            "2..101",
            "９",
        ]
        for source in invalid:
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source)
        self.assertEqual(
            parse_location(
                str(SAFE_INTEGER),
                context=SequenceContext(SAFE_INTEGER, "linear"),
            ),
            Point(SAFE_INTEGER),
        )

    def test_archived_uncertain_point_ranges_are_preserved(self) -> None:
        within = self.parse("12.21")
        self.assertEqual(within, Within(12, 21))
        self.assertEqual(
            normalize_segments(within, context=LINEAR),
            (NormalizedSegment(None, 11, 21, 1, "uncertain-point"),),
        )
        for source in ("12.12", "21.12", "<12.21", "12.>21"):
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source)

    def test_rejects_wraps_and_fuzz_errors(self) -> None:
        invalid = [
            "90..10",
            ">1..10",
            "1..<10",
            "1..>10<",
            "<1",
            ">1",
            "<1^2",
            "1^>2",
            "",
        ]
        for source in invalid:
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source, context=CIRCULAR)

    def test_ascii_spaces_are_ignored_except_inside_components_and_operator_heads(self) -> None:
        self.assertEqual(self.parse(" 1 .. 2 "), Span(Endpoint(1), Endpoint(2)))
        self.assertEqual(
            self.parse("join( 1 .. 2 , complement( 9 .. 10 ) )"),
            Join(
                (
                    Span(Endpoint(1), Endpoint(2)),
                    Complement(Span(Endpoint(9), Endpoint(10))),
                )
            ),
        )
        self.assertEqual(
            self.parse("NC_000001.11 : 3 .. 8"),
            Remote("NC_000001.11", Span(Endpoint(3), Endpoint(8))),
        )
        for source in ("1 2", "NC_000 001.11:1", "join (1,2)", "1\t..2"):
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source)

    def test_rejects_double_complement_compound_nesting_and_operator_mixing(self) -> None:
        invalid = [
            "complement(complement(1..2))",
            "join(join(1,2),3)",
            "join(order(1,2),3)",
            "order(join(1,2),3)",
            "join(complement(join(1,2)),3)",
            "join(1)",
            "order(1)",
        ]
        for source in invalid:
            with self.subTest(source=source):
                with self.assertRaises(INSDCLocationError):
                    self.parse(source)

    def test_depth_node_and_child_budgets_fail_closed(self) -> None:
        cases = [
            ("complement(1)", LocationLimits(max_depth=1), "depth"),
            ("join(1,2)", LocationLimits(max_nodes=2), "AST nodes"),
            (
                "join(1,2)",
                LocationLimits(max_compound_children=1),
                "children",
            ),
            ("join(1,2)", LocationLimits(max_children=1), "children"),
            ("123", LocationLimits(max_bytes=2), "bytes"),
        ]
        for source, limits, message in cases:
            with self.subTest(source=source, limits=limits):
                with self.assertRaisesRegex(INSDCLocationError, message):
                    self.parse(source, limits=limits)

    def test_context_and_manual_ast_construction_are_closed(self) -> None:
        for length, topology in ((0, "linear"), (True, "linear"), (1, "unknown")):
            with self.subTest(length=length, topology=topology):
                with self.assertRaises(INSDCLocationError):
                    SequenceContext(length, topology)
        with self.assertRaisesRegex(INSDCLocationError, "double complement"):
            Complement(Complement(Point(1)))
        with self.assertRaisesRegex(INSDCLocationError, "nesting"):
            Join((Join((Point(1), Point(2))), Point(3)))
        with self.assertRaisesRegex(INSDCLocationError, "uppercase accession"):
            Remote("alias", Point(1))
        with self.assertRaisesRegex(INSDCLocationError, "closed INSDC"):
            normalize_segments(object(), context=LINEAR)  # type: ignore[arg-type]


if __name__ == "__main__":
    unittest.main()
