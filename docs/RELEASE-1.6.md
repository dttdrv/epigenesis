# Epigenesis 1.6

Version 1.6.0 adds two reproducible neuroscience workflows and rewrites the
README around running and inspecting them.

## New workflows

- Neural evidence reconstructs 109 supplied CHE-1 depletion observations from
  Traets et al., Figure 6F. It links measurements to workbook cells, reconstructs
  intended donor edits, compares representations, and exports SVG, CSV and JSON.
  The offline verifier regenerates every output byte from the pinned inputs.
- Neural patterning connects an Nkx2.2 enhancer substitution and measured GLI
  binding profiles to a published four-gene regulatory model. It includes
  signal-history experiments, rescue and phenocopy, passive-reporter and
  endogenous-edit modes, and numerical refinement checks. The binding-profile
  to intracellular-affinity transfer remains an explicit biological assumption.
- Both viewers use paper-style figures: thin black axes, serif labels, no grid
  lines, and distinct point or line styles. Patterning plots retain every sampled
  value; labeled segments show contiguous sampled expression domains. SVG exports
  include the selected mode, model time, scientific qualification and provenance.
- Printing preserves the selected experiment and representation, opens scientific
  disclosures for the report, and restores their previous state afterward.

The source archive includes both examples, viewers, tests, scientific inputs
and attribution. The wheel continues to install only the compiler and its
validators. The demos require Python 3.11+ and no additional packages or network.

## Compatibility

The universal DNA translation compiler retains its common typed source boundary,
profile-native source structure and caller-supplied interpretation. Package
version changes to 1.6.0; compiler producer identities and serialized contracts
remain unchanged. Version-1 external closures remain supported alongside the
version-2 executable, digest-pinned frontend and validator path.

## Research status at the 1.6.0 release

The measured-data workflow is usable independently of predictor development.
Current neuronal DNA-edit candidates have not passed both predictive acceptance
requirements. The reserved confirmation remains unopened. A frozen training-only
pooled DNA-shape diagnostic failed its progression screens and justified no
additional fit; the ledger remains 196 fits. The
[roadmap](NEURAL_PREDICTOR_ROADMAP.md) records the result and its scope.

Natural-genome-to-brain development remains the research objective. The new
patterning experiment tests a conditional regulatory mechanism; its biological
transfer needs independent intervention evidence.

## Research and presentation update, 18 October 2026

The fixed follow-up studies are complete, with the bundled predictor unchanged.
The sequence-model research ledger records 445
fits. The two-cycle mixture has validation utility 0.07158510 versus k5's
0.03510890, and MSE 0.02862935 versus k5's 0.02862511. It has not passed both
requirements; confirmation remains unopened. The corrected quadratic-readout
receipt preserves its exact ranking tie and the original rounding discrepancy.

A new [measurement example](../examples/neural-measurement/README.md) makes the
completed RNA-repeat pilot portable. It includes all twelve original count-pair
inputs, moment archives and results, independent scoring functions, exact
arithmetic tests, source identities and methods. The command verifies relative
hashes and reproduces the saved calculations without fitting or sampling.
BIRD mean CRPS is 0.05490724 versus 0.06219211 for pooling, an 11.7% reduction
on four fixed TRAIN components. The 99% interval covers simulation error on this
cohort. A fresh DNA-normalized comparison reduces CRPS by 15.8% versus pooling
and passes both directional checks; a post hoc zero-effect forecast scores
better. The treatment-response comparison remains inconclusive.

Temporal forecasting now includes both Timestral-5M and Timestral-209M. A simple
two-feature ridge control outperforms both on the pooled-activity task. The
subsequent raw-count test holds out infection cultures and sequence families;
its four-control acceptance rule is not met. All 15 evaluation groups complete,
with 11,696 scored curves and 164 retained predictions lacking complete targets.
These are separate targets and comparisons, not a shared percentage-error scale.

The new [forecast archive](../examples/neural-forecast/README.md) includes the
frozen protocol, unchanged analysis source, count totals, coefficients, seals,
support ledgers and results. Its portable independent replay reproduces the
observable and all scores without fitting. The two temporal studies add 420
and 945 small ridge fits to their own ledgers; sequence weights remain fixed.

The README leads with the working neuroscience use case and links the complete
research history. A [presentation guide](NEUROSCIENCE_PRESENTATION.md) provides
an offline demonstration, a development map, the measured result and a proposed
controlled experiment. The [prior-art comparison](PRIOR_ART.md) now includes
neural self-construction, Gli regulation, BIRD and COMBINE archives.

Package version remains 1.6.0. Compiler contracts and predictor weights are
unchanged. Both replays are included in the source archive and CI; their NumPy
dependency does not enter the compiler wheel.

## Verification of the prepared update

Checked on macOS, 18 October 2026:

- Main suite: 504 tests, 501 passed and 3 skipped. The skips still require the
  absent sibling runtime (two) or Linux address-space limits (one).
- All 15 release checks passed, including the installed wheel and the source
  archive tested outside the checkout. Five predictor contract tests and the
  CPU inference command passed with the pinned dependencies.
- Measurement replay reproduced all twelve archived calculations; its largest
  absolute discrepancy was 3.40e-13 across the compared values. Four independent
  arithmetic tests passed. Relocation passed; altered counts, a missing moment
  archive, a reversed gain and a deliberately wrong CRPS calculation were rejected.
- Raw-count forecasting passed 23 synthetic checks before count conversion;
  12 deliberately introduced faults were detected. Independent replay verified
  all 15 evaluation groups and 30 sealed prediction files within 1.78e-15.
  Three portable-archive tests cover relocation, altered counts and a deliberately
  wrong feature calculation. The failed scientific acceptance result is retained.
- Neural evidence generation and byte replay passed. Patterning passed all
  twelve mechanism checks. Browser inspection exercised a measurement's source
  record, passive-reporter mode, model-time selection and intervention disclosure.

These are local checks. Remote CI and untested platforms require their own runs.

## Verification at the 1.6.0 release

Verified locally on macOS with Python 3.11.15 and Java 21:

- Unit suite: 502 tests, 499 passed and 3 skipped. Two skips require the absent
  sibling runtime checkout; one requires Linux address-space limits.
- All 15 release checks passed: source causality, universal translation, minimal
  example, mechanism controls, neural fate, regulation, construction, coding
  consequences, binding, boundary, contract attacks, resource/path safety,
  installed real-DNA development, wheel and source archive.
- Neural evidence generation and byte-for-byte replay passed. The patterning
  demo passed all 12 mechanism and numerical checks.
- Both viewers were inspected at desktop and 390-pixel mobile widths. Source
  inspection, comparison controls, filtering, reporter mode and playback were
  exercised. Mobile plots scroll within their containers.
- The actual patterning SVG export handler passed six independent fixture cases
  with non-integer model times, units and mode-specific legends. Reintroducing
  the old sample-index label caused the regression check to fail.

- The actual plotting script preserved values, gene order and exact sampled
  domain endpoints across both modes and nonuniform test positions. Five
  deliberately broken variants failed those checks. Print lifecycle checks also
  caught five deliberate defects, including lost disclosure state and model time.
- Both final four-page PDF reports were checked for figure layout and retained
  scientific context. Safari's Save PDF workflow also produced a readable report.

Browser file-download completion remains unverified: the in-app browser displayed
the export link, but its download-event wait timed out. Safari requested download
permission; that prompt was cancelled. SVG payload checks verify generated
content, not browser delivery. Cross-platform CI results are not implied by these
local checks.
