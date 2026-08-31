"""Bounded parser for the INSDC feature-table location language.

This module implements the closed location subset specified by INSDC Feature
Table Definition 11.4.  It preserves location syntax as an immutable AST and
performs no accession lookup, alias resolution, sequence retrieval, or
biological inference.  Callers must provide sequence bounds explicitly.

Coordinates in the AST remain 1-based and closed.  :func:`normalize_segments`
is the only conversion to 0-based half-open coordinates.  Unknown ordinary
remote coordinates remain explicitly unresolved; circular-origin syntax
requires an explicit remote length and topology.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Mapping, TypeAlias

from ._canonical import ContractError, SAFE_INTEGER


MAX_LOCATION_BYTES = 1024 * 1024
MAX_LOCATION_DEPTH = 32
MAX_LOCATION_NODES = 300_005
MAX_LOCATION_CHILDREN = 100_000
MAX_COMPOUND_CHILDREN = 100_000
MAX_REFERENCE_CONTEXTS = 100_000
MAX_ACCESSION_BYTES = 64

_REFERENCE_RE = re.compile(r"([A-Z][A-Z0-9_]*[0-9])\.([1-9][0-9]*)\Z")
_REMOTE_PREFIX_RE = re.compile(r"([A-Z][A-Z0-9_]*[0-9])\.([1-9][0-9]*):")
_SPLIT_COMPONENT_RE = re.compile(r"[A-Za-z0-9_*'\-] +[A-Za-z0-9_*'\-]")
_SPLIT_OPERATOR_RE = re.compile(r"(?:complement|join|order) +\(")


class INSDCLocationError(ContractError):
    """An INSDC location or its explicit sequence context is invalid."""


def _fail(code: str, detail: str) -> INSDCLocationError:
    return INSDCLocationError(f"{code}: {detail}")


def _coordinate(value: object, label: str) -> int:
    if type(value) is not int or not 1 <= value <= SAFE_INTEGER:
        raise _fail(
            "LOC004",
            f"{label} must be an integer in [1, {SAFE_INTEGER}]",
        )
    return value


def _validate_reference(value: object) -> str:
    if type(value) is not str or len(value) > MAX_ACCESSION_BYTES:
        raise _fail(
            "LOC003",
            f"remote reference must be at most {MAX_ACCESSION_BYTES} ASCII bytes",
        )
    match = _REFERENCE_RE.fullmatch(value)
    if match is None:
        raise _fail(
            "LOC003",
            "remote reference must be an uppercase accession with an exact version",
        )
    version = match.group(2)
    if _decimal_value(version, "remote accession version") > SAFE_INTEGER:
        raise _fail("LOC003", "remote accession version exceeds the I-JSON safe range")
    return value


def _decimal_value(token: str, label: str) -> int:
    """Convert a nonempty ASCII decimal without Python's large-int parser."""

    value = 0
    for character in token:
        digit = ord(character) - ord("0")
        if not 0 <= digit <= 9:
            raise _fail("LOC003", f"{label} is not an ASCII decimal")
        if value > (SAFE_INTEGER - digit) // 10:
            raise _fail("LOC004", f"{label} exceeds the I-JSON safe range")
        value = value * 10 + digit
    return value


@dataclass(frozen=True, slots=True)
class SequenceContext:
    """Explicit length and topology for one sequence reference."""

    length: int
    topology: str

    def __post_init__(self) -> None:
        _coordinate(self.length, "sequence length")
        if type(self.topology) is not str or self.topology not in ("linear", "circular"):
            raise _fail("LOC005", "sequence topology must be 'linear' or 'circular'")


@dataclass(frozen=True, slots=True)
class LocationLimits:
    """Per-parse ceilings, each no greater than the production hard ceiling."""

    max_bytes: int = MAX_LOCATION_BYTES
    max_depth: int = MAX_LOCATION_DEPTH
    max_nodes: int = MAX_LOCATION_NODES
    max_children: int = MAX_LOCATION_CHILDREN
    max_compound_children: int = MAX_COMPOUND_CHILDREN

    def __post_init__(self) -> None:
        ceilings = {
            "max_bytes": MAX_LOCATION_BYTES,
            "max_depth": MAX_LOCATION_DEPTH,
            "max_nodes": MAX_LOCATION_NODES,
            "max_children": MAX_LOCATION_CHILDREN,
            "max_compound_children": MAX_COMPOUND_CHILDREN,
        }
        for name, ceiling in ceilings.items():
            value = getattr(self, name)
            if type(value) is not int or not 1 <= value <= ceiling:
                raise _fail("LOC007", f"{name} must be an integer in [1, {ceiling}]")


class LocationNode:
    """Marker base for the exact closed location AST."""

    __slots__ = ()


@dataclass(frozen=True, slots=True)
class Endpoint:
    """A 1-based span endpoint and its literal INSDC fuzz marker."""

    position: int
    fuzz: str | None = None

    def __post_init__(self) -> None:
        _coordinate(self.position, "endpoint position")
        if self.fuzz is not None and (
            type(self.fuzz) is not str or self.fuzz not in ("<", ">")
        ):
            raise _fail("LOC003", "endpoint fuzz must be '<', '>', or absent")


@dataclass(frozen=True, slots=True)
class Point(LocationNode):
    position: int

    def __post_init__(self) -> None:
        _coordinate(self.position, "point coordinate")


@dataclass(frozen=True, slots=True)
class Span(LocationNode):
    start: Endpoint
    end: Endpoint

    def __post_init__(self) -> None:
        if type(self.start) is not Endpoint or type(self.end) is not Endpoint:
            raise _fail("LOC006", "span endpoints must be Endpoint values")
        if self.start.fuzz not in (None, "<") or self.end.fuzz not in (None, ">"):
            raise _fail(
                "LOC003",
                "span fuzz is only valid as '<' on the start and '>' on the end",
            )
        if self.start.position > self.end.position:
            raise _fail(
                "LOC004",
                "span start exceeds end; circular origin crossing must be explicit",
            )


@dataclass(frozen=True, slots=True)
class Between(LocationNode):
    left: int
    right: int

    def __post_init__(self) -> None:
        _coordinate(self.left, "between left coordinate")
        _coordinate(self.right, "between right coordinate")


@dataclass(frozen=True, slots=True)
class Within(LocationNode):
    """Archived single-base uncertainty over an inclusive coordinate range."""

    left: int
    right: int

    def __post_init__(self) -> None:
        _coordinate(self.left, "within left coordinate")
        _coordinate(self.right, "within right coordinate")
        if self.left >= self.right:
            raise _fail("LOC004", "within-range coordinates must increase")


AtomicLocation: TypeAlias = Point | Span | Between | Within


@dataclass(frozen=True, slots=True)
class Remote(LocationNode):
    reference: str
    location: AtomicLocation

    def __post_init__(self) -> None:
        _validate_reference(self.reference)
        if type(self.location) not in {Point, Span, Between, Within}:
            raise _fail("LOC006", "a remote reference must qualify one atomic location")


@dataclass(frozen=True, slots=True)
class Complement(LocationNode):
    location: LocationNode

    def __post_init__(self) -> None:
        if type(self.location) is Complement:
            raise _fail("LOC006", "double complement is not in the closed grammar")
        if type(self.location) not in {Point, Span, Between, Within, Remote, Join, Order}:
            raise _fail("LOC006", "complement contains an invalid location node")


def _compound_child(node: object) -> bool:
    if type(node) in {Join, Order}:
        return True
    return type(node) is Complement and type(node.location) in {Join, Order}


def _validate_compound(children: object, operator: str) -> None:
    if type(children) is not tuple:
        raise _fail("LOC006", f"{operator} children must be an immutable tuple")
    if not 2 <= len(children) <= MAX_COMPOUND_CHILDREN:
        raise _fail(
            "LOC006",
            f"{operator} requires 2..{MAX_COMPOUND_CHILDREN} children",
        )
    allowed = {Point, Span, Between, Within, Remote, Complement}
    for child in children:
        if type(child) not in allowed or _compound_child(child):
            raise _fail(
                "LOC006",
                "compound nesting and join/order mixing are not in the closed grammar",
            )


@dataclass(frozen=True, slots=True)
class Join(LocationNode):
    children: tuple[LocationNode, ...]

    def __post_init__(self) -> None:
        _validate_compound(self.children, "join")


@dataclass(frozen=True, slots=True)
class Order(LocationNode):
    children: tuple[LocationNode, ...]

    def __post_init__(self) -> None:
        _validate_compound(self.children, "order")


Location: TypeAlias = Point | Span | Between | Within | Remote | Complement | Join | Order


@dataclass(frozen=True, slots=True)
class NormalizedSegment:
    """One source-expression segment in 0-based half-open coordinates."""

    reference: str | None
    start: int
    end: int
    orientation: int
    kind: str
    start_fuzz: str | None = None
    end_fuzz: str | None = None
    bounds_status: str = "verified"

    def __post_init__(self) -> None:
        if self.reference is not None:
            _validate_reference(self.reference)
        if (
            type(self.start) is not int
            or type(self.end) is not int
            or not 0 <= self.start <= self.end <= SAFE_INTEGER
        ):
            raise _fail("LOC004", "normalized segment bounds are invalid")
        if type(self.orientation) is not int or self.orientation not in (-1, 1):
            raise _fail("LOC006", "normalized orientation must be -1 or 1")
        if type(self.kind) is not str or self.kind not in ("interval", "between", "uncertain-point"):
            raise _fail("LOC006", "normalized segment kind is invalid")
        if self.start_fuzz not in (None, "<") or self.end_fuzz not in (None, ">"):
            raise _fail("LOC006", "normalized segment fuzz is invalid")
        if self.kind == "interval" and self.start == self.end:
            raise _fail("LOC006", "interval segments must contain at least one base")
        if self.kind == "between" and (
            self.start != self.end
            or self.start_fuzz is not None
            or self.end_fuzz is not None
        ):
            raise _fail("LOC006", "between segments must be exact zero-width boundaries")
        if self.kind == "uncertain-point" and (
            self.start >= self.end
            or self.start_fuzz is not None
            or self.end_fuzz is not None
        ):
            raise _fail("LOC006", "uncertain-point segment bounds are invalid")
        if self.bounds_status not in {"verified", "unresolved"}:
            raise _fail("LOC006", "normalized segment bounds status is invalid")


def _contexts(
    references: Mapping[str, SequenceContext] | None,
) -> dict[str, SequenceContext]:
    if references is None:
        return {}
    if not isinstance(references, Mapping):
        raise _fail("LOC001", "reference contexts must be a mapping")
    if len(references) > MAX_REFERENCE_CONTEXTS:
        raise _fail(
            "LOC007",
            f"reference contexts exceed {MAX_REFERENCE_CONTEXTS} entries",
        )
    result: dict[str, SequenceContext] = {}
    for reference, context in references.items():
        reference = _validate_reference(reference)
        if type(context) is not SequenceContext:
            raise _fail("LOC001", f"context for {reference!r} must be SequenceContext")
        result[reference] = context
    return result


def _bound_atomic(
    location: AtomicLocation,
    context: SequenceContext | None,
) -> None:
    if type(location) is Point:
        coordinates = (location.position,)
    elif type(location) is Span:
        coordinates = (location.start.position, location.end.position)
    elif type(location) is Between:
        coordinates = (location.left, location.right)
    elif type(location) is Within:
        coordinates = (location.left, location.right)
    else:
        raise _fail("LOC006", "value is not an atomic location")
    if context is not None and any(value > context.length for value in coordinates):
        raise _fail("LOC005", "location coordinate exceeds its explicit sequence length")
    if type(location) is Between:
        adjacent = location.left < SAFE_INTEGER and location.right == location.left + 1
        circular_origin = (
            location.right == 1
            and context is not None
            and context.topology == "circular"
            and location.left == context.length
        )
        if not adjacent and not circular_origin:
            raise _fail(
                "LOC004",
                "between-base coordinates must be adjacent or explicit circular length^1",
            )


class _Parser:
    def __init__(
        self,
        source: str,
        context: SequenceContext,
        references: dict[str, SequenceContext],
        limits: LocationLimits,
    ) -> None:
        self.source = source
        self.context = context
        self.references = references
        self.limits = limits
        self.position = 0
        self.nodes = 0
        self.children = 0

    def parse(self) -> Location:
        location = self._location(1, inside_compound=False)
        if self.position != len(self.source):
            raise _fail(
                "LOC003",
                f"unexpected token at byte {self.position + 1}",
            )
        return location

    def _reserve_node(self) -> None:
        self.nodes += 1
        if self.nodes > self.limits.max_nodes:
            raise _fail("LOC007", f"location exceeds {self.limits.max_nodes} AST nodes")

    def _depth(self, depth: int) -> None:
        if depth > self.limits.max_depth:
            raise _fail("LOC007", f"location exceeds depth {self.limits.max_depth}")

    def _location(self, depth: int, *, inside_compound: bool) -> Location:
        self._depth(depth)
        if self.source.startswith("complement(", self.position):
            self._reserve_node()
            self.position += len("complement(")
            child = self._location(depth + 1, inside_compound=inside_compound)
            self._literal(")")
            if type(child) is Complement:
                raise _fail("LOC006", "double complement is not in the closed grammar")
            if inside_compound and type(child) in {Join, Order}:
                raise _fail("LOC006", "compound locations cannot be nested")
            return Complement(child)
        if self.source.startswith("join(", self.position):
            return self._compound(depth, "join", inside_compound)
        if self.source.startswith("order(", self.position):
            return self._compound(depth, "order", inside_compound)
        if self.position < len(self.source) and "A" <= self.source[self.position] <= "Z":
            return self._remote(depth)
        return self._atomic(self.context)

    def _compound(self, depth: int, operator: str, inside_compound: bool) -> Location:
        if inside_compound:
            raise _fail(
                "LOC006",
                "compound nesting and join/order mixing are not in the closed grammar",
            )
        self._reserve_node()
        self.position += len(operator) + 1
        children: list[LocationNode] = []
        while True:
            if len(children) >= self.limits.max_compound_children:
                raise _fail(
                    "LOC007",
                    f"{operator} exceeds {self.limits.max_compound_children} children",
                )
            self.children += 1
            if self.children > self.limits.max_children:
                raise _fail(
                    "LOC007",
                    f"location exceeds {self.limits.max_children} compound children",
                )
            child = self._location(depth + 1, inside_compound=True)
            children.append(child)
            if self._take(","):
                continue
            self._literal(")")
            break
        if len(children) < 2:
            raise _fail("LOC006", f"{operator} requires at least two children")
        immutable = tuple(children)
        return Join(immutable) if operator == "join" else Order(immutable)

    def _remote(self, depth: int) -> Remote:
        self._depth(depth)
        self._reserve_node()
        match = _REMOTE_PREFIX_RE.match(self.source, self.position)
        if match is None:
            raise _fail(
                "LOC003",
                "remote location requires an uppercase accession and exact version",
            )
        reference = f"{match.group(1)}.{match.group(2)}"
        _validate_reference(reference)
        self.position = match.end()
        context = self.references.get(reference)
        location = self._atomic(context)
        return Remote(reference, location)

    def _atomic(self, context: SequenceContext | None) -> AtomicLocation:
        self._reserve_node()
        left_fuzz = self._fuzz()
        left = self._integer("coordinate")
        if self._take(".."):
            right_fuzz = self._fuzz()
            right = self._integer("coordinate")
            span = Span(Endpoint(left, left_fuzz), Endpoint(right, right_fuzz))
            _bound_atomic(span, context)
            return span
        if self._take("^"):
            if left_fuzz is not None:
                raise _fail("LOC003", "between-base coordinates cannot be fuzzy")
            right_fuzz = self._fuzz()
            if right_fuzz is not None:
                raise _fail("LOC003", "between-base coordinates cannot be fuzzy")
            right = self._integer("between coordinate")
            between = Between(left, right)
            _bound_atomic(between, context)
            return between
        if self._take("."):
            if left_fuzz is not None:
                raise _fail("LOC003", "archived uncertain ranges cannot be fuzzy")
            right = self._integer("within coordinate")
            within = Within(left, right)
            _bound_atomic(within, context)
            return within
        if left_fuzz is not None:
            raise _fail("LOC003", "fuzz markers are only accepted on span endpoints")
        point = Point(left)
        _bound_atomic(point, context)
        return point

    def _fuzz(self) -> str | None:
        if self._take("<"):
            return "<"
        if self._take(">"):
            return ">"
        return None

    def _integer(self, label: str) -> int:
        start = self.position
        if start >= len(self.source) or not self.source[start].isdigit():
            raise _fail("LOC003", f"expected {label} at byte {start + 1}")
        if self.source[start] == "0":
            raise _fail("LOC004", f"{label} must be positive with no leading zero")
        value = 0
        while self.position < len(self.source) and self.source[self.position].isdigit():
            digit = ord(self.source[self.position]) - ord("0")
            if value > (SAFE_INTEGER - digit) // 10:
                raise _fail("LOC004", f"{label} exceeds the I-JSON safe range")
            value = value * 10 + digit
            self.position += 1
        return value

    def _take(self, token: str) -> bool:
        if self.source.startswith(token, self.position):
            self.position += len(token)
            return True
        return False

    def _literal(self, token: str) -> None:
        if not self._take(token):
            raise _fail("LOC003", f"expected {token!r} at byte {self.position + 1}")


def parse_location(
    source: str,
    *,
    context: SequenceContext,
    references: Mapping[str, SequenceContext] | None = None,
    limits: LocationLimits | None = None,
) -> Location:
    """Parse one complete INSDC 11.4 location into the closed immutable AST.

    ``references`` is an explicit accession.version-to-context mapping.  A
    remote location may be parsed without an entry, but circular-origin syntax
    and normalization require the relevant entry because neither topology nor
    length is inferred.
    """

    if type(source) is not str:
        raise _fail("LOC001", "location source must be text")
    if type(context) is not SequenceContext:
        raise _fail("LOC001", "local context must be SequenceContext")
    if limits is None:
        limits = LocationLimits()
    if type(limits) is not LocationLimits:
        raise _fail("LOC001", "limits must be LocationLimits")
    if not source:
        raise _fail("LOC003", "location source is empty")
    if len(source) > limits.max_bytes:
        raise _fail("LOC007", f"location exceeds {limits.max_bytes} bytes")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as failure:
        raise _fail("LOC003", "location syntax must be ASCII") from failure
    if any(character.isspace() and character != " " for character in source):
        raise _fail("LOC003", "location syntax may contain ASCII spaces only")
    if _SPLIT_COMPONENT_RE.search(source) is not None:
        raise _fail("LOC003", "spaces cannot split a component name or coordinate")
    if _SPLIT_OPERATOR_RE.search(source) is not None:
        raise _fail("LOC003", "an operator cannot be separated from its opening parenthesis")
    compact = source.replace(" ", "")
    if not compact:
        raise _fail("LOC003", "location source contains no descriptor")
    return _Parser(compact, context, _contexts(references), limits).parse()


def normalize_segments(
    location: Location,
    *,
    context: SequenceContext,
    references: Mapping[str, SequenceContext] | None = None,
) -> tuple[NormalizedSegment, ...]:
    """Resolve a parsed location to deterministic 0-based half-open segments.

    Segment order follows expression evaluation.  Consequently complementing
    a compound reverses child order as well as orientation.  Unknown ordinary
    remote coordinates are marked unresolved without retrieval or inference;
    remote circular-origin syntax requires explicit bounds and topology.
    """

    if type(context) is not SequenceContext:
        raise _fail("LOC001", "local context must be SequenceContext")
    reference_contexts = _contexts(references)

    def visit(
        node: LocationNode,
        active_reference: str | None,
        active_context: SequenceContext | None,
    ) -> tuple[NormalizedSegment, ...]:
        if type(node) is Remote:
            remote_context = reference_contexts.get(node.reference)
            return visit(node.location, node.reference, remote_context)
        if type(node) is Point:
            _bound_atomic(node, active_context)
            return (
                NormalizedSegment(
                    active_reference,
                    node.position - 1,
                    node.position,
                    1,
                    "interval",
                    bounds_status="verified" if active_context is not None else "unresolved",
                ),
            )
        if type(node) is Span:
            _bound_atomic(node, active_context)
            return (
                NormalizedSegment(
                    active_reference,
                    node.start.position - 1,
                    node.end.position,
                    1,
                    "interval",
                    node.start.fuzz,
                    node.end.fuzz,
                    "verified" if active_context is not None else "unresolved",
                ),
            )
        if type(node) is Between:
            _bound_atomic(node, active_context)
            boundary = (
                0
                if node.right == 1
                and (
                    active_context is None
                    or node.left == active_context.length
                )
                else node.left
            )
            return (
                NormalizedSegment(
                    active_reference,
                    boundary,
                    boundary,
                    1,
                    "between",
                    bounds_status="verified" if active_context is not None else "unresolved",
                ),
            )
        if type(node) is Within:
            _bound_atomic(node, active_context)
            return (
                NormalizedSegment(
                    active_reference,
                    node.left - 1,
                    node.right,
                    1,
                    "uncertain-point",
                    bounds_status="verified" if active_context is not None else "unresolved",
                ),
            )
        if type(node) is Complement:
            children = visit(node.location, active_reference, active_context)
            return tuple(
                NormalizedSegment(
                    segment.reference,
                    segment.start,
                    segment.end,
                    -segment.orientation,
                    segment.kind,
                    segment.start_fuzz,
                    segment.end_fuzz,
                    segment.bounds_status,
                )
                for segment in reversed(children)
            )
        if type(node) in {Join, Order}:
            result: list[NormalizedSegment] = []
            for child in node.children:
                result.extend(visit(child, active_reference, active_context))
            return tuple(result)
        raise _fail("LOC006", "value is not a closed INSDC location AST node")

    if type(location) not in {Point, Span, Between, Within, Remote, Complement, Join, Order}:
        raise _fail("LOC006", "value is not a closed INSDC location AST node")
    return visit(location, None, context)


__all__ = [
    "Between",
    "Complement",
    "Endpoint",
    "INSDCLocationError",
    "Join",
    "Location",
    "LocationLimits",
    "LocationNode",
    "MAX_ACCESSION_BYTES",
    "MAX_COMPOUND_CHILDREN",
    "MAX_LOCATION_BYTES",
    "MAX_LOCATION_CHILDREN",
    "MAX_LOCATION_DEPTH",
    "MAX_LOCATION_NODES",
    "MAX_REFERENCE_CONTEXTS",
    "NormalizedSegment",
    "Order",
    "Point",
    "Remote",
    "SequenceContext",
    "Span",
    "Within",
    "normalize_segments",
    "parse_location",
]
