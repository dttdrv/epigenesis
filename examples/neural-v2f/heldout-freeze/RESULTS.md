# Frozen Held-Out Test-Split Evaluation Results: V2F Radar on GSE188264

**Status:** Single-shot scored after cryptographic freeze verification (`FREEZE.json` SHA-256 digests verified).
**Biological acceptance:** `false` (`biological_acceptance: false`).

## Frozen SHA-256 Digests Verified Before Scoring

| File | SHA-256 |
|---|---|
| `PLAN.md` | `75ad95ba7f5e418d26ecb90ceb2afcd3e2eb52fd18cf4ff4eed886c56594e9ba` |
| `predict.py` | `acef11433c3bd86cd75c1bf9ee877e5dbb77c6aee1a2f3464ccaad03be9fc0cc` |
| `score.py` | `db2dae79409fd97352cfaa91ae699049b1f3f5c8ee4e17d52dc537623dcdcf41` |
| `predictions.json` | `298c430045403824fd316d5f081d4f5fff9acd48b16b8c61095a7819f21e1541` |

## Cohort

- **Unseen held-out `test` loci (curated 5 loci excluded):** `118` loci (`118 × 513 = 60,534` exhaustive SNVs scanned)
- **Unseen measured non-WT alleles in `test`:** `1,885` sequence perturbations (CRE/targeted/scrambled MNVs & deletions)

## Check Outcomes (Reported Verbatim, Including Null/Negative Results)

### Check H1 — Unseen Locus Cryptic Rate & Null Control

- `8 / 118` unseen test loci (`6.78%`) harbor at least one strict regime `cryptic_bifurcation` SNV (`24` cryptic SNVs total across the `118` loci; `24 / 24` are `doubly_silent` under Kreimer 0–72h `|kreimer_mpra_mean| < 0.15`).
- Top unseen cryptic locus: `chr3:20371482-20371653` (`14` cryptic SNVs).
- Null control (`n = 100` Altschul–Erickson dinucleotide shuffles + column-scrambled GLI2 PWM):
  - `dinucleotide_null_mean = 0.0` (`max = 0`, `empirical_p_ge = 0.009901`)
  - `scrambled_pwm_cryptic = 0`

### Check H2 — Passive-vs-Measured Reporter Concordance Across 1,885 Unseen Alleles

**Negative biological concordance on this GSE188264 cohort:**

- `Pearson(|max Δreporter_GLI2|, |measured α_mean|) = 0.012258`
- `Spearman(|max Δreporter_GLI2|, |measured α_mean|) = 0.019435`
- `Pearson(kreimer_mpra_mean, measured α_mean) = 0.131192`
- `Spearman(|kreimer_mpra_mean|, |measured α_mean|) = -0.138904`

The single-TF GLI2 passive reporter delta alone does **not** explain empirical variant activity across broad GSE188264 perturbations (most perturbations target other TFs or core enhancer regions).

### Check H3 — Measured Allele Outcomes by Predicted Regime Category

| Predicted Category | Count | Mean `\|measured α_mean\|` | Median `\|measured α_mean\|` | Empirically Silent (`\|α_mean\| < 0.15`) |
|---|---:|---:|---:|---:|
| `cryptic_bifurcation` | 7 | `1.258827` | `1.257264` | `0.000000` (`0 / 7`) |
| `concordant_switch` | 15 | `2.471281` | `2.260744` | `0.000000` (`0 / 15`) |
| `canalized_buffered` | 0 | `null` | `null` | `null` |
| `neutral` (under single GLI2 PWM) | 1,863 | `1.969485` | `1.346116` | `0.002147` (`4 / 1,863`) |

**Honest conclusion:** While predicted `cryptic_bifurcation` alleles show lower mean empirical MPRA effect than predicted `concordant_switch` alleles (`1.259` vs `2.471`), **none (`0 / 7`) of the measured cryptic CRE alleles are empirically silent below `0.15`** in GSE188264 (designed CRE perturbations disrupt multi-base regulatory cores and yield strong MPRA changes `~1.19–1.32`). Notice also `scram_seq57` (`siraj_score = 0.55975`, `kreimer_mpra_mean = 0.124252`, `sequence_silent = true`) illustrates why Kreimer 0–72h is the active sequence gate while Siraj (`r = 0.049`) is informational only.

All claims remain structural and model-conditional (`biological_acceptance = false`).
