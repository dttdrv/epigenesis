# Post-score zero-effect reference

The frozen transfer result is complete. It shows positive MSE contrasts against
fitted controls but uncertain correlation and panel-selection gains. Before
presenting its MSE improvement, compute the zero-effect reference on every
original study and both original aggregates. This is a post-score diagnostic,
not a new acceptance test or a change to the frozen result.

Use precisely the same held rows, equal-family weights and declared 2,000
paired family draws. Compute MSE for predictions fixed at zero, comparator
MSE minus Huber MSE and relative point reduction. Do not compute correlation
of a constant prediction. Report all six studies and both groups regardless
of sign. No fitted coefficients, scale selection, policy selection or extra
subgroups. Preserve the original scoring/RESULT.json byte-for-byte.

The independent score replay may calculate these values in a separate
ZERO_REFERENCE.json. They diagnose whether lower error than fitted baselines
also beats making no-effect predictions; they do not establish a noise floor
or fresh confirmation. No refits, new sources or sealed data.

Completed result SHA256: 360f922037d0f0cddcd9aeaead9881b5b6c5413faca7d06a6cbcea1a71df94e4
