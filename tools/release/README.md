# Epigenesis

Research platform for genome-driven neural development: genome + starting
conditions → development → neural structure → activity and learning.

**Status.** The compiler and the experiments below run and are independently
replayable. The decisive goal — predicting the neural effect of an unseen
natural DNA change better than simpler alternatives on independent data — is
**not yet met**. Every predictive test so far (six neural reporter studies,
sequence-model and motif-context variants) found at best a very small effect;
see [what is and isn't established](docs/STATUS.md).

## Quick start

Needs Python 3.11+. The two demos use only the standard library and bundled
data: no install, GPU or network.

```sh
python3 epigenesis.py evidence     # 109 published CHE-1 measurements traced to DNA edits
python3 epigenesis.py patterning   # Nkx2.2 enhancer edit in a developmental model
```

Each writes a report under `./epigenesis-out/` and prints the HTML file to open.

Score your own reporter DNA edits (270-base reference/edited pairs; format in
`examples/neural-predictor/edits.json`):

```sh
python3 -m pip install -r examples/neural-predictor/requirements.txt   # numpy, torch
python3 epigenesis.py predict --input my-edits.json --output predictions.json
```

Check this copy: `python3 epigenesis.py verify`.

## Where things are

| | |
|---|---|
| `epigenesis.py` | the commands above (`python3 epigenesis.py` lists them) |
| [`docs/STATUS.md`](docs/STATUS.md) | results and limits, in plain terms |
| below | the full project description, all examples and the compiler walkthrough |
| `examples/` | one folder per experiment, each with its own README |
| `brainc/`, [`SPEC.md`](SPEC.md) | the DNA → Development Module compiler (`pip install .` gives the `brainc` commands) |
| [`docs/research/NEURAL_PREDICTOR_ROADMAP.md`](docs/research/NEURAL_PREDICTOR_ROADMAP.md) | research record and order of work |

Compiler code: Apache-2.0 ([LICENSE](LICENSE)). Bundled data and third-party
runtimes keep their own licenses (see each example's NOTICE).
