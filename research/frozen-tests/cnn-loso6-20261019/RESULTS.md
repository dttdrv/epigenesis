# Pooled-effect CNN, leave-one-study-out (scored once, 2026-10-19)

Result: **no gain**. On the five studies without gse225815 the CNN (N1, trained on Siraj plus the other studies) gives r = 0.014 (5-95% interval -0.008 to 0.035) against 0.047 for the Siraj-trained head S2; N1 minus S2 = -0.032 (-0.060 to -0.007). Averaging them (Ceq6) gives 0.039, below S2 (-0.008, interval -0.023 to +0.006). All six: N1 0.011, S2 0.047. Utility at 50/200: N1 0.008/0.006, S2 0.046/0.034. Label not earned.

The two seeds agree poorly on held-out variants (correlation 0.15-0.51) and Siraj-TUNING loss never gets below 99% of the zero-prediction loss: the network finds little stable sequence-to-effect signal in this pool. Files: PLAN.md, FREEZE.json, derived/RESULT.json.
