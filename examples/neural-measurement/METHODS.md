# Archived RNA-repeat forecast methods

## Cohort and observations

All sources belong to the original Dev TRAIN partition. Eligible sources were
ordered by SHA-256 of `bird-predictive-pilot-v1\0dev:<source_row>`; the first four
distinct merged genomic components were selected. Their identifiers are in
`MANIFEST.json`. Admission used all three RNA repeats, so this pilot does not
constitute a prospectively independent biological test.

Counts are selected aggregate reference/alternate pairs from Rummel's DevLib
plasmid DNA and RNA samples GSM7056942, GSM7056943 and GSM7056944. The compressed
raw count source has 39,678,063 bytes and SHA-256
`c0c3fb22360c9c849075db061a7e442adf8419bdde6107d6bd2c5be390464ee9`.
The original admitted DNA and RNA exports have SHA-256
`b5bed19c9b3aa9b09602e254c3bbb3a2ff7a3f5374445e472b2c58913007ffc5`
and `aa8710df66a3333fb35ef836c31fc30626b72a5ba71b52c59ef9e6565cd797d8`.
Their complete intake pipeline is outside this saved-calculation package.

Each of the three RNA repeats is omitted in turn. Its total count is the
forecast exposure. Its alternate count enters only the score. The other two
RNA repeats and DNA counts determine the BIRD posterior. The replay checks that
the same three count pairs appear in the correct training/omitted roles.

## Forecasts

The BIRD author hierarchy is pinned in `sources.json`. The DNA alternate fraction
has prior `p ~ Beta(DNA_alt + 1, DNA_ref + 1)`. Other prior variables are
`c ~ Gamma(shape=1.1, scale=2000)`, `s ~ Gamma(shape=1.1, scale=1/3)` and independent
`w ~ Normal(0, 1)`. Values `c <= 2` contribute zero raw weight. One shared `w`
determines `q = logistic(logit(p) + s*w)` across the RNA repeats in each latent
state. Each repeat has its own RNA fraction drawn from
`Beta(1 + (c-2)*q, 1 + (c-2)*(1-q))`. A binomial draw gives the count.

The original implementation integrates the two training RNA fractions in the
importance weight, scaled by their binomial-kernel maxima. It accumulates first
and second raw weight moments in the omitted-count bins. Those saved sums and
squares are the input to this replay. The sampling driver and likelihood
implementation are not executed or requalified here.

Pooling uses `Beta(1 + sum(training_alt), 1 + sum(training_ref))`, followed by a
binomial count at the omitted total. The replay computes this Beta-binomial mass
with a recurrence. The tests independently compare it against exact finite
products. The two methods differ in priors, DNA use and heterogeneity.

## Score and numerical interval

For omitted total `n`, observed alternate count `y` and forecast CDF `F`,

```text
CRPS = (1/n) * sum[k=0..n-1] (F(k) - 1[k >= y])^2
gain = pooled_CRPS - BIRD_CRPS
```

Each of twelve jobs has equal weight, giving equal weight to each component.
`pair_loss` uses the equivalent expected-distance form of CRPS. Exact rational
tests compare both definitions, including endpoint observations.

The empirical Bernstein bound uses the original raw particle count, not an
effective-particle estimate. It encloses the weight denominator and each CDF
numerator, with a union allocation across the twelve jobs and their bins.
For variance estimate `v`, raw count `N`, cap `C` and
`L = log(4 * jobs * (n+1) / delta)`, the radius includes
`sqrt(2*v*L/N) + 7*C*L/(3*(N-1))`, plus explicit floating-point allowances.
CDF ratio bounds propagate to a CRPS interval before equal-weight aggregation.
The full arithmetic is in `scoring.py`.

The overall failure probability is `delta=0.01`. The fixed log-likelihood
allowance is `1e-8`. Coverage of the exact author posterior additionally assumes
that allowance holds throughout the prior. The earlier finite numerical checks
do not prove such a uniform bound. The interval describes integration uncertainty
for these selected observations, not biological population uncertainty.

## Fixed precision follow-up

The original 500,000-particle-per-job pilot had an inconclusive gain interval.
A separate fixed follow-up used 2,000,000 particles per job, with fresh PCG64
streams `[20261011, 3, source_row, omission]`. Model, cohort, comparator, score,
numerical allowance and acceptance rule were unchanged. All twelve jobs completed
within the fixed 45 CPU-second cap, on one thread without a GPU.

Only the follow-up interval determines its conclusion. The two runs' particles
and intervals were not pooled or intersected. The 99% statement is conditional
on the preceding history, not a simultaneous guarantee over repeated attempts.
The decision rule is positive aggregate lower gain bound for `bird_better`,
negative upper bound for `pooled_better`, otherwise `unresolved`.

The release preserves all twelve original INPUT, RESULT and moment files, plus
the aggregate result and attempt record. This portable replay verifies their
arithmetic and internal consistency. Source admission, sampling correctness and
biological generalization are separate claims.
