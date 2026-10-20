# Epigenesis

DNA compilation, neural-development simulations and experimental edit prediction.
The main demonstration changes a GLI binding site and compares four-gene dynamics
under a chosen signal history. The Huber-linear predictor provides a separate
assay-specific edit-ranking workflow. See the [current status](STATUS.md).

This quick start is copied to the root by `tools/make_release.py`. Run its
commands from the generated release folder.

## Quick start

Needs Python 3.11+. For the English/Bulgarian desktop interface, install with
[uv](https://docs.astral.sh/uv/) and start **Simulate**:

```sh
uv tool install --editable studio --with-editable .
epigenesis-studio
```

The two report commands below use only the standard library and bundled data:
no install, GPU or network.

```sh
python3 epigenesis.py evidence     # 109 published CHE-1 measurements traced to DNA edits
python3 epigenesis.py patterning   # Nkx2.2 enhancer edit in a developmental model
```

Each writes a report under `./epigenesis-out/` and prints the HTML file to open.

Score DNA edits with the same Huber-linear model used by Studio:

```sh
python3 -m pip install numpy
python3 epigenesis.py predict \
  --input studio/epigenesis_studio/example-edits.json \
  --output predictions.json --budget 3
```

Replace the example with a JSON array of `id`, `reference`, `sequence` and
optional `family`. Sequences need at least five uppercase ACGT bases; they can
have different lengths. A panel budget requires explicit families. The
[predictor guide](../../examples/neural-predictor/README.md) documents the
endpoint, input contract and evidence.

The archived 270-base Salomon models require an explicit `--activity` or
`--attention` flag. Their separate `--families` file selects ten experiments;
only attention accepts `--device`. Install
`examples/neural-predictor/requirements.txt` for either archived model.

Check the compiler and example suite: `python3 epigenesis.py verify`.
Studio's separate checks are documented in its [README](../../studio/README.md).

## Where things are

| | |
|---|---|
| `epigenesis.py` | the commands above (`python3 epigenesis.py` lists them) |
| [Status](STATUS.md) | results and scientific scope, copied to `docs/STATUS.md` |
| below | the full project description, all examples and the compiler walkthrough |
| `examples/` | one folder per experiment, each with its own README |
| `brainc/`, [`SPEC.md`](../../SPEC.md) | the DNA → Development Module compiler (`pip install .` gives the `brainc` commands) |
| [Roadmap](../../docs/research/NEURAL_PREDICTOR_ROADMAP.md) | current model summary followed by the dated research record |

Compiler code: Apache-2.0 ([LICENSE](../../LICENSE)). Bundled data and third-party
runtimes keep their own licenses (see each example's NOTICE).
