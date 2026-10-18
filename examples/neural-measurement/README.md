# Neuronal assay measurement forecast

A portable replay of a completed omitted-RNA-repeat pilot. It recomputes twelve
forecast scores from archived importance moments and an analytic pooled
Beta-binomial comparator. It uses no random draws, model fitting or new data.

From the repository root, with Python 3.12 and NumPy installed:

```sh
python examples/neural-measurement/replay.py
python -m unittest discover -s examples/neural-measurement -p 'test_*.py' -v
```

The original runtime was Python 3.12.13 with NumPy 2.5.1. The dependency is pinned
in `requirements.txt`; the compiler itself does not require NumPy.

| Fixed four-component pilot | Mean CRPS, lower is better |
|---|---:|
| BIRD repeat model | 0.05490723648 |
| Pooled Beta-binomial | 0.06219211366 |
| Relative loss reduction | 11.7135% |

The gain is 0.00728488, with conditional 99% simulation interval
[0.00436489, 0.00975665]. All four component means improve. Source 325, omission 2
favors pooling and remains included. The score measures RNA alternate-fraction
forecast distributions; it is not a DNA-edit error percentage.

`replay.py` checks the relative manifest, the complete twelve-job population,
omission-role consistency, raw moment constraints and saved score agreement.
`scoring.py` retains the six arithmetic functions from the independent local
checker. Tests compare the pooled forecast and CRPS against exact rational
references. The manifest records integrity relative to this checkout; it is not
external certification.

The exported result files retain their original bytes. The replay reproduces
the saved calculation, not the original sampling or source admission. The
[methods](METHODS.md) describe those stages and the numerical assumptions.
[Source identities](sources.json) preserve the original run and source citations.

The uncertainty interval covers Monte Carlo error conditional on four selected
TRAIN components. Admission used all RNA repeats. It is not a biological
population interval. The two models differ in priors, DNA use and heterogeneity;
the comparison does not isolate one mechanism. Broader effect and treatment
forecasts remain unresolved under their separate acceptance criteria.

Project code uses [Apache-2.0](../../LICENSE). Numeric count excerpts are
attributed to Rummel et al. and GEO GSE225815. No dataset-specific license was
found in the retained source records. No original count matrix, article text or
BIRD source code is redistributed. BIRD is Majoros et al.'s published model.

## Sources

- [Rummel et al., 2023](https://doi.org/10.1016/j.cell.2023.09.015), neuronal assay;
  [GSE225815](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225815).
- [Majoros et al., 2020; online 2019](https://doi.org/10.1093/bioinformatics/btz545), BIRD.
- [Maurer and Pontil, 2009, Theorem 4](https://arxiv.org/abs/0907.3740), empirical Bernstein bound.
- [Gneiting and Raftery, 2007, section 4.2](https://doi.org/10.1198/016214506000001437), CRPS.
