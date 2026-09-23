# Starplast

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/starplast-wordmark-white.svg">
  <img src="docs/assets/starplast-wordmark.svg" alt="Starplast — a Toxoplasma silhouette enclosing a constellation" width="520">
</picture>

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](#install)
[![Linux, macOS, Windows](https://img.shields.io/badge/Platforms-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)](#install)
[![Qt / PyQt6](https://img.shields.io/badge/GUI-PyQt6-41CD52?logo=qt&logoColor=white)](https://einarolafsson.github.io/starplast/guide.html)
[![MIT license](https://img.shields.io/github/license/EinarOlafsson/starplast?color=3DA639)](LICENSE)

[![Documentation on main](https://github.com/EinarOlafsson/starplast/actions/workflows/docs.yml/badge.svg?branch=main&event=push)](https://github.com/EinarOlafsson/starplast/actions/workflows/docs.yml?query=branch%3Amain)
[![Python API](https://img.shields.io/badge/API-reference-007EC6)](https://einarolafsson.github.io/starplast/API.html)
[![GitHub issues](https://img.shields.io/github/issues/EinarOlafsson/starplast?logo=github)](https://github.com/EinarOlafsson/starplast/issues)
[![spaCR integration](https://img.shields.io/badge/spaCR-integrated-7B61A8)](https://github.com/EinarOlafsson/spacr)

Starplast is a desktop app for exploring gene evidence in *Toxoplasma gondii* and
*Plasmodium falciparum*. It brings expression, fitness screens, localization,
protein features, interactions, and literature annotations into one place.

Use it to look up a gene, investigate hits from a screen, compare groups of genes,
and choose candidates for follow-up experiments. The bundled maps contain 8,140
*T. gondii* genes and 5,720 *P. falciparum* genes, viewed separately.

[User guide](docs/guide.md) · [Python API](docs/API.md) ·
[Dataset catalogue](docs/datasets.md) · [Changes](CHANGELOG.md)

## Install

Python 3.10 or newer is required. The desktop app needs a display and OpenGL.

```bash
pip install starplast
starplast
```

From a checkout:

```bash
git clone https://github.com/EinarOlafsson/starplast.git
cd starplast
pip install -e .
starplast
```

The default installation runs analyses on the CPU. For an NVIDIA GPU on Linux
x86_64 with a CUDA 12 compatible driver:

```bash
pip install "starplast[gpu]"
```

GPU dependencies are large and optional. `starplast-install-gpu` can also inspect
your driver and offer an installation command. macOS and Windows use CPU analysis;
the 3D display still uses OpenGL.

The built gene tables and graphs ship with the package. Browsing them works offline.
Protein structures are downloaded when requested and cached locally. Rebuilding
the bundled data requires the original source datasets.

## Explore a screen

1. Choose the organism under **File → Species**.
2. Search for a gene ID. For *T. gondii*, **File → Import data** adds your screen results.
3. Colour and filter the map using the measurements you want to compare.
4. Click a gene to read its evidence and inspect each type of relationship separately.
5. Select a group with a lasso or brush, then export the gene list for follow-up.

[![A rotating gene map with CDPK1 and its connections illuminated](docs/screenshots/map_rotation.gif)](docs/screenshots/map_rotation.png)

The *T. gondii* map with CDPK1 (`TGME49_301440`) selected. The selected gene and
its five neighbours emit light; lines show attention-corrected literature
co-mention links. Colours show compartments. [View a still image](docs/screenshots/map_rotation.png).

Each point is a gene. Its position comes from an embedding of selected features;
nearby points have similar inputs, but proximity alone does not demonstrate a
shared function or physical interaction. Colours, relationship edges, and the
evidence panel provide the context needed to interpret the map. Missing evidence
is shown separately from measured values.

The analysis panel lets you change feature sets, build embeddings, cluster genes,
and evaluate recovery of labels held out from the input. Search scores help
prioritize candidates; they are not experimental validation. See the
[user guide](docs/guide.md) for controls and analysis settings.

## Working with spaCR

[spaCR](https://github.com/EinarOlafsson/spacr) handles microscopy and image-based
screen analysis. Its Starplast launcher installs and opens this app in a separate
Python environment. Starplast provides a place to explore the biological evidence
around those screen results. Export a gene-level table from spaCR and import it
into Starplast; launching the app does not transfer results automatically.

## Data and reproducibility

**[Full table of included data and source links →](docs/datasets.md)**

The catalogue lists all 128 registered datasets and computed layers, with their
measurements, coverage, publication references, and links to source data or inputs.
Coverage differs by organism and assay;
absence from a literature search does not establish that a gene has never been studied.

```bash
python -m starplast.paths       # show data locations and missing cache files
python -m starplast.build_graph # rebuild from available source datasets
```

| Variable | Purpose |
|---|---|
| `STARPLAST_CACHE` | Override the directory containing built gene tables and graphs |
| `STARPLAST_DATA` | Locate raw datasets used for rebuilding |
| `STARPLAST_STATE` | Override the directory for saved runs, annotations, and downloads |

Methods are described in [MATERIALS_AND_METHODS.md](MATERIALS_AND_METHODS.md).
[HANDOFF.md](HANDOFF.md) contains the development history and earlier design decisions.
The [repository review](docs/repository-review.md) describes the current architecture
and its limitations.

## Development

```bash
pip install -e ".[dev,docs]"
QT_QPA_PLATFORM=offscreen pytest tests/test_app_smoke.py tests/test_docstrings.py -q
python scripts/build_docs.py
```

Open `docs/site/index.html` for the guide and generated API reference.
Some tests need source datasets, CUDA, or a working OpenGL context; see
[development and releases](docs/releases.md) for the release checks.

Develop on `nightly` and merge checked changes into `main`. To release, update all
package versions together with `python scripts/release.py bump 0.43.0` and update
the changelog before merging. A version increase on `main` triggers checks,
builds, PyPI publishing, and a GitHub release. Ordinary merges do not publish.
The initial PyPI account setup is described in the release guide.

The source code is distributed under the [MIT license](LICENSE).
Source datasets and third-party artwork retain their own licenses and attribution.
