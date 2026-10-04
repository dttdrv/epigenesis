# Local neural construction and activity

This example develops cells, neurites and physical contacts from one founder
and a shared **artificial** recipe. The final cells and connections are not
encoded or imported. A separate activity assay stimulates each developed cell
and removes its outgoing contacts to test transmission.

It is an intermediate research result, not a natural-genome interpreter or a
biologically equivalent digital brain. The 480-base sequence is an explicit
encoding of 15 floating-point parameters. No claim about intelligence follows
from these experiments.

## Run

Use Python 3.11+, a JDK supporting `javac --release 11` (tested locally with
Temurin 21), and the standard `patch` utility. Both `java` and `javac` must be
on `PATH`. No network request or package installation occurs during the run.

```sh
python3.11 examples/neural-construction/construction.py --output /tmp/my-construction --duration 3
python3.11 tests/verify.py neural-construction
```

The output directory must not exist. `--recipe` accepts an exact 480-byte
uppercase ACGT file; without it, the executable writes its declared default
recipe. `--seed` selects the random initial conditions. Failed runs retain
partial evidence and return a nonzero exit status.

Outputs include:

- `compiled/`: exact source, caller interpretation, target, linked Development
  Module and producer-isolated validation. The module contains one founder and
  the shared program, with no edge schema or imported morphology.
- `backend/`: extracted Java source, applied patch, compiler logs and class
  identities. Only `.java` files from the pinned archive are compiled; bundled
  upstream binaries are never executed. The execution API retains the expected
  manifest hash from this build rather than trusting a mutable manifest alone.
- `development-run/`: environment, lineage, local growth inputs/actions,
  physical subdivision/merger events, per-step geometry, contact formation,
  final reciprocal connections, runtime log and audit report.
- `activity.json`: every-cell stimulation, outgoing-contact lesion and zero
  transmission controls on the final recorded graph.
- `viewer.html`: offline playback of recorded geometry and activity trials.
- `experiment.json`: source/module/runtime/activity evidence hashes. These
  record local provenance; they are not an external certification or signature.

The CLI validates and consumes the linked program before development. The
lower-level `develop` function also accepts a `Program` directly for numerical
experiments; such a call alone does not establish a source-to-runtime chain.

## Mechanism and units

Each progenitor grows toward its division-volume threshold, divides equally,
and partitions its inherited division signal. Once the signal is at or below
the differentiation threshold, it creates one axonal and one dendritic growth
cone. The shared parameters determine speed, guidance gain, branching hazard,
taper and stopping diameter. Growth sees its local direction and guidance
concentration gradient. A Gaussian guidance field is an imposed environment,
not a modeled signaling tissue. The spatial mesh begins from random noncell
points; these are numerical substrate nodes, not imported neurons.

CX3D uses micrometers and hours here: speed is μm/h, volume rate is μm³/h,
branching is a hazard per hour, and taper is μm/h. Division signal and threshold
are artificial dimensionless quantities. Branch daughters each receive the
parent diameter divided by √2, preserving cross-sectional area by model choice.

Nonterminal processes form boutons or spines. Synapses require distinct cells,
reciprocal bouton/spine links, the actual spatial-neighbor/contact test and
attachment positions within the declared reach. The audit checks attachments
against their named cylinder surfaces. Recognition can be disabled. Its random
stream is separate from the morphology stream, so paired runs have identical
history before their first physical contact; later bonds can affect mechanics.

Formation excludes candidates whose stored axial coordinate falls outside a
cylinder. During passive shaft deformation, attachments follow a declared
uniform-strain rule: their material fraction `h/L` stays fixed. The code moves
each shared coordinate array once, including independent spring endpoints,
without changing owners or spring resting lengths. This applies when either
the distal mass or its parent moves; active tip addition keeps existing axial
coordinates, and mesh changes use their separate geometric transforms.

This rule makes the [upstream fixed-attachment intent](https://tilde.ini.uzh.ch/~amw/seco/cx3d/Tutorial.pdf)
explicit; it is not experimentally calibrated tissue mechanics. Earlier runs
left sites beyond a compressed shaft, and their failing evidence is retained.
The trace still reports `retained_sites_outside_segments`, so invalid geometry
cannot disappear behind reciprocal logical links. The activity assay uses the
retained graph and does not establish a biological attachment law.

The Python API supports a timed guidance intervention through
`cue_switch=(step, new_cue_center)`. The test checks identical prefixes, unchanged
local position/noise at intervention, and the resulting gradient/direction
change. Separate tests cover unseen recipes/seeds, translated initial geometry,
conserved division volume, resource exhaustion and rejection of forged traces.

The trace audit replays lineage, active growth modules, declared local growth
calculations, topology changes, surface attachment and summary closure. It does
not independently solve all CX3D mechanics, prove authenticity against a
coherently fabricated complete history, or establish order-independent physics.
The backend runs one synchronous simulation in a fresh bounded JVM.

## Activity assumptions

The assay implements the zero-current delta-input leaky integrate-and-fire
equations documented by [NEST](https://nest-simulator.readthedocs.io/en/stable/models/iaf_psc_delta.html):
exact exponential decay between events, simultaneous-input summation, threshold
and reset, and discarded input during the absolute refractory period. It does
not invoke NEST or claim equivalence to its entire simulator.

Rest/reset −70 mV, threshold −55 mV, membrane time constant 10 ms and refractory
period 2 ms are declared assay parameters. Each retained contact receives an
assumed excitatory 15 mV jump and 1 ms delay. These choices test graph
transmission; construction predicts neither synaptic efficacy nor excitatory
or inhibitory identity. Parallel contacts retain their multiplicity. The
expected response includes only paths that arrive within the assay's 20 ms
horizon. Unstimulated, blocked-transmission and outgoing-contact lesion
controls must agree with the declared mechanism. Analytical tests separately
check leak, pulse integration, threshold, refractory boundaries and event order.

## Source and modifications

The physical substrate is the official [CX3D 0.03 archive](https://tilde.ini.uzh.ch/~amw/seco/cx3d/cx3d-0.03.zip),
908,697 bytes, SHA-256
`41bc3580f253709e23fbdea21395f33c4e60091a761483e36aa76de6280d2c34`.
See the [2009 paper](https://doi.org/10.3389/neuro.10.025.2009) and
[`vendor/provenance.json`](vendor/provenance.json).

The explicit source patch removes unused JDK-internal imports, supports
headless display initialization, makes sphere volume/diameter conversions
consistent through `Math.PI`, restores a scoped contact RNG, observes physical
subdivision/merging, updates neurite coordinate axes after soma motion, and
rejects new contact proposals outside a cylinder's current axial extent.
It also transfers attachments and spring endpoints together during subdivision
and compatible merging, transforms shared coordinate arrays once, and corrects
the axial offset in a straight merge. Public mesh coarsening retains two
segments if their diameters differ or their local axes differ beyond one ULP at
unit-vector scale. This permits normalization roundoff while preserving finer
geometry. Tests check positions, identity, ownership, volume,
substance quantity and spring rest length, including independent springs and
asymmetric attachment populations. The former implementation could shift a
spring by 18 μm, duplicate attachments or crash during a merge. Unequal-diameter
and rotated-frame probes check that incompatible geometry stays intact. These
are direct public-operation tests: the pinned physics scheduler does not reach
its merge branch in ordinary runs. Its separate biological retraction path is
unused by this consumer and remains outside these checks. Passive compression,
stretch, soma motion and both branch daughters now have distinct motion probes,
with active growth as a counterexample to unconditional coordinate rescaling.
Passive deformation preserves substance quantity and the existing constant-
diameter volume/length relation; remeshing preserves volume. Invalid material
lengths, including unchanged collapsed geometry, fail clearly. Child-deformation
failure releases the parent's lock. The corrections retain their failing witnesses and the
uniform-strain assumption remains biologically uncalibrated.
Physics and diffusion remain enabled.

The archive predates current Java toolchains. It still emits upstream
unchecked/deprecation warnings: normal compilation reports 12 missing-annotation
warnings plus unchecked/deprecation notes; a detailed lint audit found 171
warnings in total. They are retained in logs, not suppressed. Broader upstream
maintenance and remote platform validation remain unfinished.

The complete 2013 cortical G-code program was not recovered. This is a new
local construction program using the published substrate, not a reproduction
of the [2013 cortical model](https://doi.org/10.1371/journal.pcbi.1003173).
Genetic developmental programming has prior art; this example makes no claim
that the broad idea is new.

CX3D retains its GPL-3.0-or-later license and copyright notices in the archive
and [`vendor/LICENSE.CX3D`](vendor/LICENSE.CX3D). Its modifications remain covered
by that license. New example source is Apache-2.0; combined Java distribution
must preserve the applicable GPL terms. The compiler wheel contains no Java
simulator. The source archive includes this separately licensed example.

The remaining scientific gates include a measured natural-sequence mechanism,
independent developmental and intervention data, neural-type and efficacy
validation, numerical/physical calibration, and organism/region coverage.
`biological_acceptance` remains false.
