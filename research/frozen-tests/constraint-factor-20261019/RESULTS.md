# Conservation × direction (exploratory, negative)

Idea: predicted effect = P(variant matters) × direction, with "matters" taken from evolutionary conservation (phyloP100way, hg19, mean over ±10 bp) and direction from the Siraj-trained head S2.

Outcomes had already been seen, so this was exploratory and no frozen plan was written. Result: nothing to freeze.

| study | n | r(abs effect, phyloP) | r(effect, S2) | r(effect, S2 × conservation) |
|---|---|---|---|---|
| lee | 5882 | 0.004 | 0.019 | 0.009 |
| mcafee | 1667 | -0.016 | 0.092 | 0.095 |
| matoba | 4000 | 0.014 | 0.022 | 0.032 |
| gse221384 | 745 | 0.075 | -0.014 | -0.062 |
| weiss | 9229 | 0.007 | 0.115 | 0.097 |
| mean | | 0.017 | 0.047 | 0.034 |

Conservation does not predict which variants have large reporter effects (mean r 0.017), and weighting S2 by it lowers the mean r from 0.047 to 0.034. Sign accuracy of S2 among the 20% most conserved variants is 50.7%, no better than overall. gse225815 is excluded (count noise, see the precision test).

Files: `scripts/phylop_par.py` (remote bigWig read, 13 shards), `scripts/explore.py`, `derived/phylop_*of*.npz`, `derived/explore.txt`.
