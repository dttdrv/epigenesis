# Thermodynamic transfer test

Tests the one assumption the [neural-circuit](../neural-circuit/README.md)
scanner has never checked against cell data: that binding physics computed from
in-vitro PWMs transfers to regulatory activity in living cells. A saturating,
non-negative activator/quencher model whose only sequence input is a PWM
partition sum over 532 JASPAR 2024 motifs is fitted to the Kreimer et al. (2022)
0–72 h neural-induction reporter data and scored on 118 held-out families
against the repository's frozen baselines and 20 shuffled-PWM nulls.

**Outcome: the frozen decision rule is not met.** Real PWMs beat all 20 shuffled
nulls, and the model beats zero effect, but it does not beat `k5` or the
shuffled-null average with an interval above zero, and it is worse than the free-weight motif-temporal predictor
(0.2686 against 0.2521). The same site sums used linearly do better (0.2611).
See [RESULTS.md](RESULTS.md); the [plan](PLAN.md) and [hashes](FREEZE.json) were
fixed before the model was fitted.

## Run

Python 3.11+ with `numpy` and `scipy` (`requirements.txt`). No network.

```sh
python3 examples/neural-thermo/thermo.py select --output /tmp/thermo-select
python3 examples/neural-thermo/thermo.py confirm \
  --selection examples/neural-thermo/selection.json --output /tmp/thermo-confirm
python3 -m unittest discover -s examples/neural-thermo -p 'test_thermo.py' -v
```

`select` fits the frozen grid on training families and scores validation only
(about 8 minutes). `confirm` refuses to run if the plan, code or selection
differ from `FREEZE.json`, refits on train + validation, scores the test
families once and writes the report, TF strengths and predictions (about 7
minutes). `result/` holds the run reported in RESULTS.md.

## Files

| File | Role |
|---|---|
| `PLAN.md`, `FREEZE.json` | Frozen protocol, gates, and hashes of plan, code and selection |
| `sites.py` | Strand-averaged site odds → saturating occupancy sums; PWM-shuffle null |
| `thermo.py` | Model, analytic gradient, L-BFGS-B fit, grid selection, bootstrap |
| `confirm.py` | The single test scoring, nulls, controls, gates |
| `selection.json` | Validation results for all 60 grid points and the chosen one |
| `result/` | `report.json`, `tf-strengths.csv`, `strengths.npz`, `test-predictions.csv` |
| `data/motif_names.json` | JASPAR 2024 names for the 532 bundled motifs (CC BY 4.0) |

## Scope

Retrospective: the test families were examined by earlier experiments. One
assay, reporter constructs, published preprocessing. Nothing here validates the
neural-tube circuit or any fate prediction.

Sources: [Kreimer et al. (2022)](https://doi.org/10.1038/s41467-022-28659-0),
GSE188264; [JASPAR 2024](https://doi.org/10.1093/nar/gkad1059).
