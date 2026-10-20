# Attribution

The compact Epigenesis PyTorch scorer uses the rational feed-forward function
from [Interlace](https://github.com/MisulOrg/Interlace/tree/b7aeae7a3364be0f834ac2777acb887d7d2859cc)
by Misul. Its normalization, numerical rescaling, sequence input, distance bias
and scalar output are implemented for this reporter task. Interlace is licensed
under Apache-2.0; the license text is included at the repository root.

Training effects originate from [Salomon et al., 2026](https://doi.org/10.64898/2026.07.16.738760)
and the [80K-Analysis repository](https://github.com/kircherlab/80K-Analysis/tree/a686f1d552f97ebcc0224057b5a835002f7a3dc8).
The package contains Epigenesis's fitted weights and aggregate results, not
the source measurement tables. The upstream MIT notice is retained in
[LICENSE.Salomon](LICENSE.Salomon).

The separate linear activity-transfer model uses the same study's admitted
physical reporter activities as additional TRAIN supervision. Its development
comparison explored a latent activity score followed by a nonlinear map,
motivated by [MAVE-NN](https://doi.org/10.1186/s13059-022-02661-7) and
[SQUID](https://doi.org/10.1038/s42256-024-00851-5). TRAIN selection favored the
linear form. This is an assay-specific predictive experiment, not a reproduction
of those methods or an inferred biochemical mechanism.
