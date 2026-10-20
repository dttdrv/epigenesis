# Fixed Huber transfer result

The frozen Huber head has mean correlation 0.048969 across the five primary studies, versus 0.046548 for historical S2. Its paired correlation gains over JL, K5 and S2 all have intervals crossing zero. Transfer superiority and improved panel selection remain unestablished.

The primary raw-MSE aggregate improves over each fitted comparator, including the same-training JL and K5 heads. This is an assay-scale diagnostic: target definitions and scales differ across studies. The separately declared post-score zero-effect reference gives MSE 0.09529676 versus Huber 0.10437080: Huber is 9.52% worse than zero. All five primary per-study points and gain intervals favor zero. This diagnostic is separate from the original frozen comparisons, and prevents the fitted-baseline improvement from being presented as useful absolute prediction.

No coefficients, calibration or ranking policies were fitted. Every native motif feature was recomputed before scoring and matched the old cache exactly, for all 22,217 held rows. The target/order check reproduced all twelve historical S2/S3 correlations.

## Populations and per-study correlation

| Study | Edits | Families | Huber | Joint-linear | K5 | S2 | S3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| lee | 5,916 | 2,340 | 0.015199 | 0.012672 | 0.025422 | 0.018508 | 0.013374 |
| mcafee | 1,682 | 676 | 0.093491 | 0.096379 | 0.099003 | 0.091026 | 0.094681 |
| matoba | 4,013 | 97 | 0.028585 | 0.031777 | 0.038416 | 0.022548 | 0.026299 |
| gse225815 | 611 | 349 | 0.037475 | 0.042533 | 0.050521 | 0.048818 | 0.047169 |
| gse221384 | 748 | 571 | -0.015748 | -0.015852 | -0.024127 | -0.014291 | -0.020470 |
| weiss | 9,247 | 6,930 | 0.123317 | 0.115821 | 0.091206 | 0.114947 | 0.112301 |

GSE225815 is excluded from the primary five following the pre-existing count-noise correction. It remains in the fixed six-study sensitivity analysis. Matoba has only 97 families, so its 200-edit utility is unavailable. Primary U200 uses Lee, McAfee, GSE221384 and Weiss; sensitivity U200 additionally uses GSE225815.

## All frozen aggregate contrasts

Each interval below is the 5th to 95th percentile of 2,000 paired whole-family-copy draws. The 5th percentile is the one-sided 95% lower bound. Each eligible study contributes equal weight. Positive contrasts favor Huber.

### primary_five

| Endpoint | Comparator | Gain | Lower | Upper |
|---|---|---:|---:|---:|
| r | JL | 0.00080922 | -0.00569406 | 0.00846772 |
| r | K5 | 0.00298454 | -0.01078660 | 0.01731020 |
| r | S2 | 0.00242125 | -0.00598364 | 0.01093761 |
| r | S3 | 0.00373169 | -0.00448338 | 0.01241074 |
| mse | JL | 0.00844641 | 0.00724092 | 0.00979632 |
| mse | K5 | 0.00266938 | 0.00030125 | 0.00456923 |
| mse | S2 | 0.00506613 | 0.00422601 | 0.00593487 |
| mse | S3 | 0.00653055 | 0.00548700 | 0.00767011 |
| u10 | JL | 0.01773366 | -0.02837360 | 0.04273915 |
| u10 | K5 | 0.01778479 | -0.09097124 | 0.08894985 |
| u10 | S2 | 0.00233669 | -0.03710360 | 0.05215619 |
| u10 | S3 | 0.03443382 | -0.02182559 | 0.06156790 |
| u50 | JL | -0.00276876 | -0.02293752 | 0.01306332 |
| u50 | K5 | 0.01569634 | -0.01253731 | 0.04570141 |
| u50 | S2 | 0.00439417 | -0.01859290 | 0.02399385 |
| u50 | S3 | -0.00153613 | -0.02606039 | 0.02007833 |
| u200 | JL | 0.00648951 | -0.01002937 | 0.01747271 |
| u200 | K5 | 0.01742428 | -0.00150942 | 0.03422813 |
| u200 | S2 | -0.00177904 | -0.01435595 | 0.01472527 |
| u200 | S3 | 0.00095720 | -0.01489301 | 0.01553572 |
### all_six

| Endpoint | Comparator | Gain | Lower | Upper |
|---|---|---:|---:|---:|
| r | JL | -0.00016867 | -0.00654748 | 0.00718273 |
| r | K5 | 0.00031283 | -0.01313502 | 0.01414947 |
| r | S2 | 0.00012731 | -0.00798572 | 0.00787228 |
| r | S3 | 0.00149413 | -0.00656310 | 0.00936807 |
| mse | JL | 0.00847227 | 0.00703450 | 0.00999183 |
| mse | K5 | 0.00223554 | -0.00003483 | 0.00423366 |
| mse | S2 | 0.00480312 | 0.00352841 | 0.00620209 |
| mse | S3 | 0.00604916 | 0.00478388 | 0.00734877 |
| u10 | JL | 0.01376892 | -0.04710112 | 0.03429081 |
| u10 | K5 | 0.01959145 | -0.11973322 | 0.09129279 |
| u10 | S2 | -0.01908492 | -0.07280083 | 0.03417457 |
| u10 | S3 | 0.00766269 | -0.06200626 | 0.04407170 |
| u50 | JL | -0.01104448 | -0.04037474 | 0.01069407 |
| u50 | K5 | 0.00729756 | -0.02992958 | 0.03301763 |
| u50 | S2 | -0.00789921 | -0.04792832 | 0.02977252 |
| u50 | S3 | 0.01004341 | -0.02649018 | 0.04063602 |
| u200 | JL | 0.00468690 | -0.01471117 | 0.01714934 |
| u200 | K5 | 0.01256954 | -0.01068605 | 0.02756811 |
| u200 | S2 | -0.00555141 | -0.02072081 | 0.01526178 |
| u200 | S3 | 0.00359833 | -0.01334631 | 0.02103630 |

MSE gain means comparator MSE minus Huber MSE. Correlation and utility gain mean Huber minus comparator. Utility is observed signed reporter effect per selected edit. Positive mean correlation alone does not establish superiority to a baseline.

## Boundaries and execution

H, JL and K5 use the same Siraj FIT population. Historical S2 and S3 used FIT plus TUNING; they are descriptive controls. All external outcomes were previously exposed in earlier experiments. These within-study intervals do not account for adaptive model search, shared study lineages or uncertainty across laboratories.

The inherited 2-kb exclusions are not comprehensive sequence/exposure admission. Several inputs are fetched 150-base genomic windows; Weiss uses supplied 230-base strings. Reconstructed features do not prove identity to every physical assay insert. These reporter targets are not measured neural development. No replicate halves or noise ceiling were used.

Prediction took 55.77 seconds and scoring took 2.40 seconds, on one CPU thread, with no spending, downloads or refits. Both stages completed within their fixed 180-second caps.

Original Siraj acceptance remains six of seven, with its ten-edit comparison against k5 unresolved. This diagnostic neither changes that gate nor opens Siraj RESERVE or Salomon confirmation. The predictor goal remains active.

Artifacts: [protocol](PLAN.md), [freeze](FREEZE.json), [prediction seal](prediction/SEAL.json), [complete scores](scoring/RESULT.json), and [post-score zero-reference protocol](POSTSCORE_ZERO_REFERENCE.md).

Independent replay passed: all 18 native H/JL/K5 prediction vectors and 12 historical vectors match exactly; all 30 full rankings, 351 score summaries and 156 paired contrasts were verified. Maximum score difference was 8.89e-16. The replay took 6.17 seconds without refitting. [Replay receipt](reviews/replay/REPLAY_RECEIPT.json).

Result SHA256: 360f922037d0f0cddcd9aeaead9881b5b6c5413faca7d06a6cbcea1a71df94e4

## Separate post-score zero-effect reference

| Study | Zero MSE | Huber MSE | Huber excess MSE |
|---|---:|---:|---:|
| lee | 0.06259210 | 0.07391969 | 18.10% |
| mcafee | 0.02768657 | 0.03332791 | 20.38% |
| matoba | 0.01637243 | 0.02562741 | 56.53% |
| gse225815 | 0.54197803 | 0.55026006 | 1.53% |
| gse221384 | 0.34787667 | 0.36121438 | 3.83% |
| weiss | 0.02195605 | 0.02776460 | 26.46% |

The primary zero-minus-Huber MSE gain is -0.00907403, with 5th/95th quantiles [-0.01212116, -0.00661120]. The all-six sensitivity result is also worse than zero by 5.27%. [Complete zero reference](ZERO_REFERENCE.json). No noise ceiling is inferred.
