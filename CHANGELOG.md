# Changelog

## 0.42.0

- Simplify the map to individual genes; remove Galaxy and Orthogroup summary modes.
- Add a Starplast logo and window icon.
- Rewrite the README and add a user guide and Python API guide.
- Publish generated API documentation through GitHub Pages.
- Add settings help and document application callbacks.
- Fix importing screen tables whose identifier column is already named `gene_id`.
- Fix Plasmodium gene selection and link its evidence panel to PlasmoDB.
- Make CUDA dependencies optional through `starplast[gpu]`.
- Include SVG artwork and compressed sequence tables in wheels; exclude saved embeddings.
- Add synchronized version bumps, build checks, and automatic PyPI publishing.

## 0.41.0

Existing development baseline before the packaging and documentation overhaul.
Earlier changes are recorded in [HANDOFF.md](HANDOFF.md) and the Git history.
