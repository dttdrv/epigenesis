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

## Subsequent repository update

The experimental neuronal edit predictor is now available as a separate
[Python 3.12 example](../examples/neural-predictor/README.md), with pinned PyTorch
and NumPy dependencies, frozen weights, CPU/Metal inference and aggregate
research receipts. The compiler retains its standard-library-only installation.

The latest two-cycle mixture has validation selection utility 0.07158510,
103.9% above its comparator, and MSE 0.02862935, 0.01480% above its strict limit.
The joint gate remains unmet after 280 biological fits. The new cohort's
confirmation remains unopened. A completed TRAIN-only uncertainty audit covers
21,781 edits and identifies error increases from uncalibrated corrections in
all four uncertainty groups. The README now presents current results in a
compact table; the roadmap retains the experiment history. Package version
remains 1.6.0; this update does not declare predictive acceptance.

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
