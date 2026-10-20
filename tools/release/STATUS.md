# Status

Plain-language summary of what the bundled experiments do and do not show.
Details and numbers are in the [roadmap](research/NEURAL_PREDICTOR_ROADMAP.md).

## Works and is replayable
- Compiler: DNA source → explicit interpretation → sealed Development Module, with a separate validator.
- Evidence report: 109 published measurements traced to workbook cells and intended DNA edits.
- Patterning, construction, fate, coding-consequence examples: executable mechanisms with checks.
- Measurement forecast: omitted-repeat forecasts improve on a fixed four-component pilot.

## Not established
- Predicting the downstream neural effect of an unseen natural DNA alteration, beating simpler models on independent data.
- The bundled edit predictor (activity-transfer and attention models) improves point estimates on held-out edits, but every confirmation interval includes zero.

## Cross-study tests (2026-10-19)
Models trained on one set of neural reporter studies were scored on held-out variants from others, always against simple baselines:

| Test | Result |
|---|---|
| k-mer/motif model trained on three studies, tested on the others | no signal |
| six studies, Siraj-trained head and combinations | r ≈ 0.05 (under 1% of variance); combining models adds nothing |
| Sei chromatin-profile features | no signal (r ≈ 0); blend no better |
| motif-context (squared) model | same as the simpler head |

Measurement noise caps are 0.6–0.98, so the models, not the data, are the limit.
Next: a model of a different kind, then application to a developmental question.
