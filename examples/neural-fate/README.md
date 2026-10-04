# Neural-fate mechanism experiment

This example computes a neural progenitor's Pax6, Olig2, and Nkx2.2 regulatory
dynamics from a compiled parameter recipe. It reproduces specified behaviors
of published ordinary differential equation models. It uses Python 3.11+ and
the standard library.

```sh
demo_dir=$(mktemp -d)
python3.11 examples/neural-fate/experiment.py --output "$demo_dir/neural-fate"
python3.11 tests/verify.py neural-fate-experiment
```

Run from the repository or unpacked source distribution on supported Linux or
macOS hosts. The output directory must not exist. Failure returns a nonzero
exit code and preserves any files already produced for inspection. The entire
experiment directory is not an atomic publication; its compiler bundles use
Epigenesis's atomic, no-clobber publication.

## Recipe and execution

The recipe is an explicitly artificial encoding of 16 model parameters. Each
IEEE 754 binary64 value is serialized little-endian; each byte becomes four
bases, with two-bit digits 0/1/2/3 mapped to A/C/G/T in most-significant-digit
order. The parameter order is:

```text
alpha beta gamma h1 h2 h3 h4 h5 k1 k2 k3
ncrit_p ocrit_p ncrit_o ocrit_n pcrit_n
```

There are exactly 512 bases. Production rates may be zero; all other parameters
must be positive. All values must be finite. This encoding contains no natural
gene annotation, cell list, adjacency matrix, precomputed fate, or trajectory.
The regulatory interactions are specified by the published equations in the
target consumer. They are not discovered from the recipe.

`compile_recipe` decodes the recipe as the caller's interpretation. It compiles
one progenitor with constant parameters, initial state `(alpha/k1, 0, 0)`, and
an attached update rule. The target contract pins the exact `experiment.py`
bytes. `load_compiled` runs the producer-isolated Development Module validator,
requires the supported target and operations, then compares the linked tensors
against the recipe. A correctly sealed but wrong interpretation is rejected
at this last check. The compiler itself does not infer the interpretation.

`simulate` consumes those linked parameters and initial state. A classical
fourth-order Runge-Kutta integrator splits steps at every signal change.
`dt` is a maximum step size. Nonfinite values, negative states, unrepresentable
time schedules, and trajectories above one million integration steps fail
explicitly. No clipping or replacement result hides numerical failure.

To test another valid parameter recipe from Python, call `encode_recipe`,
`compile_recipe`, and `load_compiled`, then pass the loaded values and a
`[(duration, signal), ...]` schedule to `simulate`. The CLI runs the fixed
reproduction panel below.

## Source equations and protocols

Let `P`, `O`, and `N` denote Pax6, Olig2, and Nkx2.2. `G` is the model's external
signal input. All concentrations and times below are model units.

```text
dP/dt = alpha / (1 + (N/ncrit_p)^h1 + (O/ocrit_p)^h2) - k1*P
dO/dt = beta*G / ((1+G)*(1 + (N/ncrit_o)^h3)) - k2*O
dN/dt = gamma*G / ((1+G)*(1 + (O/ocrit_n)^h4 + (P/pcrit_n)^h5)) - k3*N
```

[Balaskas et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3267043/),
equations 1-3 and [supplementary Table S2](https://researchonline.lshtm.ac.uk/id/eprint/2067511/1/mmc1.pdf),
provide the first parameter set. Signal activation has Hill exponent one.
The second pair comes from [Panovska-Griffiths et al. (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3565701/),
Appendix C and Figures 2 and 5. The Figure 5 settings use `h3=h4=h5=5` and
contrast `ncrit_p=0.9` with `1.1`. Exact parameter tables and retrieved source
digests are in [models.json](models.json). Papers are not redistributed.

| Condition | Parameters and signal schedule | Observed criterion |
|---|---|---|
| Naive versus primed | 2012; `40@G=1` versus `20@G=5, 20@G=1` | Final Nkx2.2 below versus above 1 |
| Signal withdrawal | 2012; `20@G=5, 20@G=0` | Final Nkx2.2 below `1e-7` |
| Double gene loss | 2012; `alpha=beta=0` versus baseline, `20@G=0.5` | Nkx2.2 above versus below 1 |
| Steady versus oscillating | 2013 Figure 5 pair; `500@G=5` | Late Olig2 amplitude below `1e-4` versus above `0.1` |

The signal-history schedules are chosen to demonstrate the model's mechanism.
They are not a numerical reproduction of the paper's 18-hour chick experiment.
The threshold of one follows the 2012 model's high/low expression convention.
The small residual and oscillation thresholds are numerical acceptance
criteria, not biological measurement uncertainties. No measured trajectories
are used for independent fitting or holdout prediction here.

## Evidence and falsification

`report.json` records every pass/fail criterion, endpoints, parameter-document
and implementation hashes, compiled bundle identities, protocols, sample
counts, and the complete trajectory file's SHA-256. `trajectories.csv` contains
every coarse integration step. Each model directory contains the exact recipe,
interpretation chain, compiled module, and validation report.

The oscillator must persist in both the 200-300 and 400-500 time windows and
have at least three peaks in the latter. The full 500-unit trajectory is
repeated at `dt=0.01` instead of `0.02`. Maximum state disagreement at shared
times must remain below `1e-4`; mean periods must agree within the conservative
sampling bound `2*(0.02+0.01)`. These checks distinguish persistent oscillation
from a transient or a coarse-step artifact.

The regression tests additionally use exact derivatives, independent
exponential solutions, a separately transcribed midpoint integrator at several
times, gene loss and exact recipe restoration, malformed inputs, numerical
boundaries, source tampering, coherently substituted interpretation tensors,
and no-clobber output. A successful run prints `published-model-reproduced`.
The focused verification gate prints `NEURAL-FATE-EXPERIMENT-PASSED`.

This is a deterministic, fixed-network model reproduction. It does not infer
biological meaning from natural DNA, build neurons or synapses, simulate
spiking, or validate a one-to-one brain. The additional requirements for a
biological claim are in [the acceptance protocol](../../docs/NEURAL_MECHANISM_ACCEPTANCE.md).
