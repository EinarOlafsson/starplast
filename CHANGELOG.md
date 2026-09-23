# Changelog

## 0.42.1

- Use the GitHub README as the PyPI project description.
- Use full URLs for README artwork, the rotating gene map, and documentation links.

## 0.42.0

- Simplify the map to individual genes; remove Galaxy and Orthogroup summary modes.
- Adopt the Toxoplasma constellation logo, SVG wordmarks, and window icon.
- Rewrite the README and add a user guide and Python API guide.
- Publish generated API documentation through GitHub Pages.
- Add settings help and document application callbacks.
- Fix importing screen tables whose identifier column is already named `gene_id`.
- Fix Plasmodium gene selection and link its evidence panel to PlasmoDB.
- Publish the application and bundled data as one PyPI project, `starplast`.
- Make CUDA dependencies optional through `starplast[gpu]`.
- Include SVG artwork and compressed sequence tables in wheels; exclude saved embeddings.
- Add synchronized version bumps, build checks, and automatic PyPI publishing.
- Develop on `nightly`; publish version increases on `main` to PyPI and GitHub Releases.
- Link all 128 registered datasets and computed layers to their sources.
- Show a 1440×1080, 30 fps gene-map rotation with selected-gene lighting in the README.

## 0.41.0

Existing development baseline before the packaging and documentation overhaul.
Earlier changes are recorded in [HANDOFF.md](HANDOFF.md) and the Git history.
