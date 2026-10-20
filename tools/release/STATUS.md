# Current delivery — 20 October 2026

The interactive entry point is [Epigenesis Studio](../../studio/README.md),
with English and Bulgarian interfaces. The compiler package remains version 1.6.0.

## DNA-site simulation

Simulate accepts a supported single substitution in the nine-base Nkx2.2 GLI
binding site, compiles its nucleotides into a validated Development Module,
then derives the binding-profile ratio from the measured table. A published
four-gene model produces reference, variant, reporter and reversion trajectories
under the selected signal history. Transferring the measured ratio to
intracellular affinity is an explicit biological assumption.

This is a working conditional mechanism experiment. The broader
genome → development → neural structure → activity → learning path remains
the research objective; its separate examples are not one validated biological model.

## Delivered predictor

Studio and `epigenesis.py predict` use the frozen Huber-linear checkpoint for
the Siraj SK-N-SH 24-hour reporter effect. It combines 5-mer counts with 532
compressed motif potentials, scores paired sequences and can select one edit
per family. It runs with NumPy on CPU and does not train during use.

On reused internal tuning data (16,986 edits, 1,759 families), equal-family
MSE is 0.28766713 versus k5's 0.29661090: a 3.02% reduction. Six of seven
original requirements pass; the ten-edit gain lower bound is −0.150.
On five external neural reporter studies, mean correlation is 0.0490 and
aggregate MSE is 9.52% above zero-effect prediction. Those results support an
experimental assay-specific tool; transferable prediction remains unestablished.

The earlier activity-transfer and attention models remain available explicitly
for reproducing the separate Salomon assay results. Their confirmation is spent;
Siraj RESERVE remains sealed.

## Evidence and reproduction

- [Predictor guide](../../examples/neural-predictor/README.md): current commands, input contract and archived models.
- [Technical paper](../../docs/paper/epigenesis-technical-paper.pdf) ([Bulgarian](../../docs/paper/epigenesis-technical-paper.bg.pdf)): mechanisms, results and scientific qualifications.
- [Frozen research](../../research/README.md): thirteen cross-study and variant-prioritization tests, with plans and hash verification.
- [Roadmap](../../docs/research/NEURAL_PREDICTOR_ROADMAP.md): current delivery followed by the dated experiment history.

The evidence reconstruction, developmental mechanisms and measurement forecasts
retain their own endpoints and controls. Reported standard errors and repeat
agreement are uncertainty diagnostics; they do not establish a common
irreducible-noise floor for these assays.
