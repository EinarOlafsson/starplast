# Development and releases

## Local checks

Use a virtual environment and install the development and documentation tools:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,docs]"
python scripts/release.py check
QT_QPA_PLATFORM=offscreen pytest tests/test_release.py tests/test_packaging.py tests/test_readme.py tests/test_docstrings.py tests/test_paths.py tests/test_app_smoke.py tests/test_app_polish.py tests/test_discover_cli.py tests/test_importer.py -q
python scripts/build_docs.py
```

On Windows, activate `.venv\Scripts\activate` instead. On Linux, Qt may need
`libegl1`, `libgl1`, `libdbus-1-3`, `libxkbcommon-x11-0`, and `libxcb-cursor0` from
the system package manager. `QT_QPA_PLATFORM=offscreen` is for tests, not normal use.

The full test suite also exercises source-dataset ingestion, rendering, GPU
backends, and long searches. Those tests need the corresponding source data or
hardware. The release checks above cover the installable package and entry points;
they are not a validation of all scientific results.

## Package layout

| Distribution | Contents |
|---|---|
| `starplast-core` | Python modules, console commands, built gene data, and icons |
| `starplast` | Installs the exact matching core version; forwards `gpu` and `ingest` extras |
| `starplast-gpu` | Compatibility alias for the matching `starplast[gpu]` version |

The core name remains stable because spaCR checks its installed metadata. CUDA is
optional for all entry points. The core wheel includes cached data for offline
browsing and excludes local saved embeddings. Source datasets are not packaged.

## One-time PyPI setup

Sign in at [PyPI publishing](https://pypi.org/manage/account/publishing/) and create
a pending Trusted Publisher for each of `starplast-core`, `starplast`, and
`starplast-gpu` with these values:

| Field | Value |
|---|---|
| Owner | `EinarOlafsson` |
| Repository | `starplast` |
| Workflow filename | `release.yml` |
| Environment | `pypi` |

For an existing project, add the same publisher in that project's Publishing
settings. See [PyPI's setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).
The workflow exchanges its GitHub identity for short-lived upload credentials;
no permanent PyPI token is required in GitHub secrets.

Create the GitHub environment `pypi` and restrict deployments to the `main` branch.
To publish automatically, leave required reviewers unset. Under repository
**Settings → Pages**, choose **GitHub Actions** as the source. The documentation
workflow builds pull requests and deploys pushes to `main`.

## Publish a version

```bash
python scripts/release.py bump 0.43.0
python scripts/release.py check
```

Use a version greater than the current one. The bump command updates the core,
both metapackages, exact internal dependency pins, and `starplast.__version__`.
Update `CHANGELOG.md`, run the release checks, commit, and push to `main`.

The release workflow compares the version with the repository state before the
push. An unchanged version does not publish; a downgrade fails. A version increase
runs the checks, builds wheels and source distributions, verifies their contents,
and publishes all three distributions. Tags are not required.

For the first upload or to retry a partially completed upload, manually run
**Publish Python packages** on `main`. Existing files are skipped; PyPI does not
allow an uploaded version to be replaced. Use a new version for changed code.
If publication failed before any upload, rerunning the original job preserves the
same tested artifacts.

## Inspect release artifacts

```bash
python -m build
python -m build packaging/starplast --outdir dist
python -m build packaging/starplast-gpu --outdir dist
python scripts/check_wheel.py dist
python -m twine check --strict dist/*
```

Start with an empty `dist/` directory so files from older versions are not included.
The wheel check requires both organism caches, the application icon and diagrams,
and the compressed sequence tables. It rejects saved embeddings, unconditional
CUDA dependencies, and files above PyPI's default upload limit.

The desktop installers in `packaging/` are separate from PyPI. Their existing
workflow runs manually or on a `v*` tag.
