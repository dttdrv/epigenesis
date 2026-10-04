# Acceptance of a DNA-directed neural mechanism

The research claim to test is: **a specified developmental recipe predicts the
effect of an unseen intervention through a declared biological mechanism.**
Acceptance applies to the named organism, tissue, developmental stage, input
range, and measured observables. It cannot establish a one-to-one digital brain.

Epigenesis currently validates and compiles a caller-supplied interpretation.
An authentic source, correctly sealed module, or changed artifact hash does not
establish that DNA caused a predicted phenotype. This is the distinction the
executable controls below enforce.

## Run the controls

```sh
python3.11 tests/verify.py neural-mechanism-controls
```

Success prints `NEURAL-MECHANISM-CONTROLS-PASSED` followed by an explicit statement
that biological development is not accepted by this gate. The gate runs in CI
and is included in ordinary unittest discovery. It adds no production runtime.

### What is executed

The positive control is an **artificial mathematical encoding** in A/C/G/T:
the first two bases encode an amplitude `A`, the next two an affinity `K`, and
the last eight are neutral by construction. A pair is a base-4 integer plus
one, with A/C/G/T denoting digits 0/1/2/3. There is no claim that natural DNA
uses this encoding or that these bases are biologically neutral.

For a nonnegative external signal `S`, the stipulated equilibrium observable is
`y = A*S/(K+S)`. The oracle uses exact rational arithmetic on the experimental
parameters; it does not use the candidate's decoder. The positive candidate
decodes the sequence and evaluates the algebraically equivalent expression
`A/(1+K/S)`, with zero at `S=0`. This checks a declared contract, not a biological
theory or dynamical simulation. All quantities are in artificial units.

Each panel has four generated parameter pairs, three signal levels, and eight
conditions per pair/level: baseline, increased K, exact allele rescue, mechanism
rescue, mechanism phenocopy, and three neutral substitutions. For mechanism
rescue, the mutated sequence stays fixed while K is clamped to its original
value. For phenocopy, the original sequence stays fixed while K is clamped to
the mutant value. These interventions test the declared mediator beyond simple
sequence replay.

The positive control runs all 96 conditions, repeats them in reverse execution
order, and runs another 96-condition panel with new generated recipes. Candidates
receive only the sequence, signal, and optional K intervention. Case names,
expected values, paths, and FASTA labels are not candidate inputs. Temporary
paths and record labels vary independently of condition. This is an interface
restriction, not a security sandbox against malicious Python code.

Every finite response passes through real source compilation, request/response
binding, typed lowering, on-disk bundle publication, and producer-isolated
validation. The test reloads the published module and decodes the tensor
referenced by its initializer. It compares that value to the oracle with
`1e-12` relative and absolute arithmetic tolerances. These tolerances are not
experimental uncertainty or biological effect thresholds.

### Required rejection controls

| Candidate or failure | Discriminating observation |
|---|---|
| Frozen caller interpretation | Mutation changes sequence identity but leaves a wrong observable; both artifacts can still validate |
| Hash-derived phenotype | A changed hash cannot substitute for quantitative agreement |
| Correct direction, 10% wrong magnitude | Direction alone fails the expected response curve |
| Environment-blind response | Zero and nonzero signal conditions require different responses |
| Neutral-region-sensitive response | Ordinary recipes pass; a neutral substitution must trigger rejection |
| Intervention-blind response | Natural alleles pass; mediation or phenocopy must trigger rejection |
| Lookup of the first panel | Previously stored observations fail on the second panel |
| NaN or infinity | A nonfinite prediction cannot pass a numerical comparison |

The public, deterministic panels are regression controls. They are not blinded
biological holdouts, and a program written to reproduce the formula can pass.
Black-box endpoint tests cannot distinguish that program from a local growth
mechanism. This gate validates the assay's discrimination and compiler transport;
it does not validate neuronal development, spiking activity, or biological DNA
interpretation.

## Executable published mechanism

The [neural-fate example](../examples/neural-fate/README.md) adds a caller-owned
parameter interpreter and a target consumer. It compiles an explicitly
artificial 512-base recipe into model parameters and initial state, validates
the module, and integrates published Pax6/Olig2/Nkx2.2 equations. No final
phenotypes or trajectories are stored in the recipe.

```sh
python3.11 tests/verify.py neural-fate-experiment
```

The example reproduces signal-history dependence, double gene loss, and
steady-versus-oscillating dynamics. Analytic limits, a separate midpoint
integrator, and timestep refinement check numerical behavior. A coherently
changed caller interpretation passes the compiler's structural checks but
fails the consumer's recipe-to-parameter comparison. That distinction is
intentional: valid provenance cannot establish a correct interpretation.

These are published-model reproductions, with investigator-selected schedules
in model units. They do not establish a natural sequence-to-mechanism mapping,
independent biological prediction, cell growth, or circuit construction.

## Measured sequence-effect experiment

The [neural-regulation example](../examples/neural-regulation/README.md) learns
reporter perturbation effects from actual enhancer sequences and published
neural-induction measurements (Kreimer et al. 2022, GSE188264). It merges
overlapping or sequence-identical families before splitting, selects on
validation families, retains zero estimates, and evaluates frozen test families
against zero-effect and dinucleotide baselines. Its sequence provider and
consumer compile and replay the fitted predictions outside `brainc/`.

```sh
python3.11 tests/verify.py neural-regulation
```

This command validates the software. The separate first real-data benchmark
**failed predictive acceptance**: its improvement over the dinucleotide model
was uncertain. Time-specific prediction also lacked demonstrated benefit.
The study's reporter estimates were pooled and normalized before our split, so
these are retrospective family holdouts rather than independent preparation
holdouts. The example establishes an executable measured-sequence experiment;
it does not meet the developmental mechanism requirements below or calibrate
the published ODE's production parameters.

## Scientific acceptance protocol

All rows below are required for a claim about a biological developmental
mechanism. A missing row means **not established**, even when software tests
pass. These are evidence requirements, not implemented or prechecked results.

| Requirement | Evidence and falsification test |
|---|---|
| Define the claim before fitting | Specify species, stage, tissue, initial state, environmental inputs, observables, units, measurement times, parameter bounds, and exclusions. Freeze the model and analysis before revealing test outcomes. |
| Establish sequence-to-mechanism mapping | Identify actual loci and variants; cite measured or independently validated effects on expression, binding, channels, or other modeled processes. Audit every translation step. Separately label measured, fitted, inferred, and assumed parameters. A nucleotide-to-arbitrary-number table cannot satisfy this row. |
| Withhold interventions, not adjacent samples | Separate training and test by biological preparation and intervention family; reserve variants and signal histories. Track donor/batch dependence. Do not treat neurons from one preparation as independent replicates or refit after a failed test. |
| Measure development over time | Record the developmental variables actually claimed, at prespecified times: for example cell fate and local signals for a fate mechanism, or lineage, position, growth, and connectivity for circuit construction. State exclusions explicitly. Check transitions as well as the final readout. A final graph alone cannot establish a construction mechanism. |
| Audit local construction | Declare allowed cell inputs and interactions. Record enough state to replay updates from initial conditions. Inspect source and runtime access for imported adjacency, per-cell answer tables, or hidden target-state feedback. Changing unused distant state must not immediately change a local update when the specified interaction rules forbid it. |
| Test causal interventions | Use baseline, causal perturbation, independently supported neutral/sham controls, and experimental rescue. Intervene at the proposed mediator while holding the variant fixed; induce the predicted phenotype through that mediator with the original allele. Score predicted intermediate states and downstream observations. Do not assume exact biological rescue. |
| Test environmental dependence | Withhold doses, timing, durations, and initial states within the declared domain. Include cases where a mutation should have little effect and cases where it should matter. Compare full response trajectories or distributions, not only signs or selected pictures. |
| Establish numerical reliability | Test time-step and spatial-resolution refinement, conservation/bounds where the equations require them, deterministic replay, and ensembles under declared random seeds. Preserve failures. Numerical error must be smaller than the prespecified distinguishable biological effect. |
| Compare competing explanations | Evaluate a static reconstruction, a simpler fitted response model, and mechanism ablations under the same data access and fitting budget. Select interventions where their predictions differ. Equally good predictions mean the proposed mechanism remains unresolved. |
| Quantify predictive accuracy | Prespecify primary endpoints and acceptable discrepancies using measurement uncertainty and a scientifically meaningful effect. Use preparation-level replication and uncertainty intervals; handle multiple endpoints and missing data in the frozen analysis. Report all held-out interventions, calibration, residuals, and comparisons. No universal percentage or arbitrary p-value makes a brain model valid. |

For stochastic models, compare distributions and prespecified contrasts across
independent biological preparations. Paired simulation seeds can reduce Monte
Carlo noise in an intervention contrast, but they do not create new biological
replicates. Choose the ensemble and experiment size from the target precision;
do not stop at the first favorable result.

A compact first biological target is a neural progenitor fate mechanism with
measured signal-history and transcription-factor perturbations. For example,
Balaskas et al. study the Pax6/Olig2/Nkx2.2 regulatory network that interprets
Sonic Hedgehog signaling. Reproducing a supplied regulatory network would test
that mechanism; it would not establish inference of the network from raw DNA.
[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3267043/)

Gene-like recipes for local neural self-construction have prior art. Zubler et
al. describe local developmental primitives and a simulated network grown from
a progenitor. The useful research question is which additional interventions
our model can predict, with which evidence and limits.
[Primary paper](https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2011.00057/full)

## Evidence record for a biological result

Preserve exact sequence and dataset identities, licenses/provenance, train/test
assignment, executable and parameter identities, initial conditions, seed list,
environment schedules, interventions, complete trajectories, measurement
extraction code, numerical convergence results, baseline/ablation predictions,
and the frozen scoring rule. Report compiler validation separately from
mechanistic and predictive results. An independently auditable record is needed
to distinguish a successful prediction from a reconstructed answer.
