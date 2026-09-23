# Starplast

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/assets/starplast-wordmark-white.svg">
  <img src="https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/assets/starplast-wordmark.svg" alt="Starplast — a Toxoplasma silhouette enclosing a constellation" width="520">
</picture>

[![Platforms](https://img.shields.io/badge/Platforms-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)](#install)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](#install)
[![Qt / PyQt6](https://img.shields.io/badge/GUI-PyQt6-41CD52?logo=qt&logoColor=white)](https://einarolafsson.github.io/starplast/guide.html)
[![Latest release](https://img.shields.io/github/v/release/EinarOlafsson/starplast?label=Release&logo=github)](https://github.com/EinarOlafsson/starplast/releases/latest)
[![GitHub issues](https://img.shields.io/github/issues/EinarOlafsson/starplast?logo=github)](https://github.com/EinarOlafsson/starplast/issues)
[![GitHub source](https://img.shields.io/badge/GitHub-Source-181717?logo=github)](https://github.com/EinarOlafsson/starplast)

[![PyPI version](https://img.shields.io/pypi/v/starplast?logo=pypi&logoColor=white)](https://pypi.org/project/starplast/)
[![Total PyPI downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fsql-clickhouse.clickhouse.com%2F%3Fuser%3Ddemo%26query%3DSELECT%2Bif%2528count%2528%2529%2B%253D%2B0%252C%2B%2527awaiting%2Bdata%2527%252C%2BtoString%2528sum%2528count%2529%2529%2529%2BAS%2Btotal_text%252C%2Bif%2528count%2528%2529%2B%253D%2B0%252C%2B%2527awaiting%2Bdata%2527%252C%2BtoString%2528sumIf%2528count%252C%2Bdate%2B%253E%253D%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2B-%2B30%2BAND%2Bdate%2B%253C%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2529%2529%2529%2BAS%2Bmonth_text%2BFROM%2Bpypi.pypi_downloads_per_day%2BWHERE%2Bproject%2B%253D%2B%2527starplast%2527%2BFORMAT%2BJSON&query=%24.data%5B0%5D.total_text&label=downloads&color=brightgreen&cacheSeconds=86400)](https://clickpy.clickhouse.com/dashboard/starplast)
[![PyPI downloads over 30 complete days](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fsql-clickhouse.clickhouse.com%2F%3Fuser%3Ddemo%26query%3DSELECT%2Bif%2528count%2528%2529%2B%253D%2B0%252C%2B%2527awaiting%2Bdata%2527%252C%2BtoString%2528sum%2528count%2529%2529%2529%2BAS%2Btotal_text%252C%2Bif%2528count%2528%2529%2B%253D%2B0%252C%2B%2527awaiting%2Bdata%2527%252C%2BtoString%2528sumIf%2528count%252C%2Bdate%2B%253E%253D%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2B-%2B30%2BAND%2Bdate%2B%253C%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2529%2529%2529%2BAS%2Bmonth_text%2BFROM%2Bpypi.pypi_downloads_per_day%2BWHERE%2Bproject%2B%253D%2B%2527starplast%2527%2BFORMAT%2BJSON&query=%24.data%5B0%5D.month_text&label=downloads+%2830d%29&color=brightgreen&cacheSeconds=86400)](https://clickpy.clickhouse.com/dashboard/starplast)
[![PyPI download rank over 30 days](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fsql-clickhouse.clickhouse.com%2F%3Fuser%3Ddemo%26param_package_name%3Dstarplast%26param_days%3D30%26query%3DWITH%2B%2528%2BSELECT%2Bsum%2528count%2529%2BFROM%2Bpypi.pypi_downloads_per_day%2BWHERE%2Bproject%2B%253D%2B%257Bpackage_name%253AString%257D%2BAND%2Bdate%2B%253E%253D%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2B-%2B%257Bdays%253AUInt16%257D%2BAND%2Bdate%2B%253C%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2B%2529%2BAS%2Bdownloads%2BSELECT%2Bdownloads%2BAS%2Bpackage_downloads%252C%2BcountIf%2528n%2B%253E%253D%2Bdownloads%2529%2BAS%2Brank%252C%2Bcount%2528%2529%2BAS%2Btotal_packages%252C%2B100.0%2B%252A%2Brank%2B%252F%2BnullIf%2528total_packages%252C%2B0%2529%2BAS%2Bpercentile%252C%2Bif%2528%2Btotal_packages%2B%253D%2B0%2BOR%2Bdownloads%2B%253D%2B0%252C%2B%2527no%2Bdata%2527%252C%2Bconcat%2528%2B%2527top%2B%2527%252C%2BtoString%2528ceil%25281000.0%2B%252A%2Brank%2B%252F%2BnullIf%2528total_packages%252C%2B0%2529%2529%2B%252F%2B10%2529%252C%2B%2527%2525%2527%2B%2529%2B%2529%2BAS%2Bmessage%2BFROM%2B%2528%2BSELECT%2Bproject%252C%2Bsum%2528count%2529%2BAS%2Bn%2BFROM%2Bpypi.pypi_downloads_per_day%2BWHERE%2Bdate%2B%253E%253D%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2B-%2B%257Bdays%253AUInt16%257D%2BAND%2Bdate%2B%253C%2BtoDate%2528now%2528%2527UTC%2527%2529%2529%2BGROUP%2BBY%2Bproject%2B%2529%2BFORMAT%2BJSON&query=%24.data%5B0%5D.message&label=PyPI+rank+%2830d%29&color=brightgreen&cacheSeconds=86400)](https://clickpy.clickhouse.com/dashboard/starplast)

[![Tests on main](https://github.com/EinarOlafsson/starplast/actions/workflows/checks.yml/badge.svg?branch=main&event=push)](https://github.com/EinarOlafsson/starplast/actions/workflows/checks.yml?query=branch%3Amain)
[![Documentation on main](https://github.com/EinarOlafsson/starplast/actions/workflows/docs.yml/badge.svg?branch=main&event=push)](https://github.com/EinarOlafsson/starplast/actions/workflows/docs.yml?query=branch%3Amain)
[![User guide](https://img.shields.io/badge/Guide-Using%20Starplast-4A9EFF)](https://einarolafsson.github.io/starplast/guide.html)
[![Python API](https://img.shields.io/badge/API-reference-007EC6)](https://einarolafsson.github.io/starplast/API.html)
[![Cite Starplast](https://img.shields.io/badge/Cite-CITATION.cff-8A2BE2)](https://github.com/EinarOlafsson/starplast/blob/main/CITATION.cff)
[![MIT license](https://img.shields.io/github/license/EinarOlafsson/starplast?color=3DA639)](https://github.com/EinarOlafsson/starplast/blob/main/LICENSE)
[![spaCR integration](https://img.shields.io/badge/spaCR-integrated-7B61A8)](https://github.com/EinarOlafsson/spacr)

Starplast is a desktop app for exploring gene evidence in *Toxoplasma gondii* and
*Plasmodium falciparum*. It brings expression, fitness screens, localization,
protein features, interactions, and literature annotations into one place.

Use it to look up a gene, investigate hits from a screen, compare groups of genes,
and choose candidates for follow-up experiments. The bundled maps contain 8,140
*T. gondii* genes and 5,720 *P. falciparum* genes, viewed separately.

[User guide](https://github.com/EinarOlafsson/starplast/blob/main/docs/guide.md) · [Python API](https://github.com/EinarOlafsson/starplast/blob/main/docs/API.md) ·
[Dataset catalogue](https://github.com/EinarOlafsson/starplast/blob/main/docs/datasets.md) · [Changes](https://github.com/EinarOlafsson/starplast/blob/main/CHANGELOG.md)

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

[![A rotating gene map with CDPK1 and its connections illuminated](https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/screenshots/map_rotation.gif)](https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/screenshots/map_rotation.png)

The *T. gondii* map with CDPK1 (`TGME49_301440`) selected. The selected gene and
its five neighbours emit light; lines show attention-corrected literature
co-mention links. Colours show compartments. This animation uses the pre-0.43
layout; newly opened maps use the recorded, balanced feature recipe.
[View a still image](https://raw.githubusercontent.com/EinarOlafsson/starplast/main/docs/screenshots/map_rotation.png).

Each point is a gene. Its position comes from an embedding of selected features;
nearby points have similar inputs, but proximity alone does not demonstrate a
shared function or physical interaction. Colours, relationship edges, and the
evidence panel provide the context needed to interpret the map. Missing evidence
is shown separately from measured values.

The analysis panel lets you change feature sets, build embeddings, cluster genes,
and evaluate recovery of labels held out from the input. Search scores help
prioritize candidates; they are not experimental validation. See the
[user guide](https://github.com/EinarOlafsson/starplast/blob/main/docs/guide.md) for controls and analysis settings.

## Explore, predict and compare

The **Tools** menu offers three starting points:

- **Explore a gene** searches identifiers, symbols and descriptions, then shows the source evidence.
- **Predict a trait** compares held-out predictions before making calls for unlabelled genes. Related genes stay in the same evaluation group; reports retain missing values, calibration status and exact input settings.
- **Compare a screen** joins a gene-level CSV/TSV to the existing evidence and keeps unresolved identifiers visible.

The package includes frozen protein sequence representations for 8,064 Toxoplasma
proteins and AF3 summary features for 1,210 exactly mapped genes. Sequence features
are optional inputs to prediction; local structure files can also be indexed for
the structure viewer.

In the [0.43 benchmark](https://github.com/EinarOlafsson/starplast/blob/main/docs/benchmark-0.43.md),
boosted trees with sequence and structure features reached 56.2% localization
accuracy across held-out orthogroups. UMAP neighbours performed less well, and
several rare classes remained poorly recovered. These are model hypotheses,
not experimentally established annotations.

[Guided workflows](https://github.com/EinarOlafsson/starplast/blob/main/docs/workflows.md) ·
[Local AF3 structures](https://github.com/EinarOlafsson/starplast/blob/main/docs/structures.md) ·
[Protein sequence features](https://github.com/EinarOlafsson/starplast/blob/main/docs/protein-sequences.md)

## Working with spaCR

[spaCR](https://github.com/EinarOlafsson/spacr) handles microscopy and image-based
screen analysis. Its Starplast launcher installs and opens this app in a separate
Python environment. Starplast provides a place to explore the biological evidence
around those screen results. Export a gene-level table from spaCR and import it
into Starplast; launching the app does not transfer results automatically.

## Data and reproducibility

**[Full table of included data and source links →](https://github.com/EinarOlafsson/starplast/blob/main/docs/datasets.md)**

The catalogue lists all 130 registered datasets and computed layers, with their
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

Methods are described in [MATERIALS_AND_METHODS.md](https://github.com/EinarOlafsson/starplast/blob/main/MATERIALS_AND_METHODS.md).
[HANDOFF.md](https://github.com/EinarOlafsson/starplast/blob/main/HANDOFF.md) contains the development history and earlier design decisions.
The [repository review](https://github.com/EinarOlafsson/starplast/blob/main/docs/repository-review.md) describes the current architecture
and its limitations.

## Development

```bash
pip install -e ".[dev,docs]"
QT_QPA_PLATFORM=offscreen pytest tests/test_app_smoke.py tests/test_docstrings.py -q
python scripts/build_docs.py
```

Open `docs/site/index.html` for the guide and generated API reference.
Some tests need source datasets, CUDA, or a working OpenGL context; see
[development and releases](https://github.com/EinarOlafsson/starplast/blob/main/docs/releases.md) for the release checks.

Develop on `nightly` and merge checked changes into `main`. To release, update all
package and runtime versions together with `python scripts/release.py bump 0.43.0` and update
the changelog before merging. A version increase on `main` triggers checks,
builds, PyPI publishing, and a GitHub release. Ordinary merges do not publish.
The initial PyPI account setup is described in the release guide.

The source code is distributed under the [MIT license](https://github.com/EinarOlafsson/starplast/blob/main/LICENSE).
Source datasets and third-party artwork retain their own licenses and attribution.
