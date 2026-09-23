# Development and releases

## Branches and commits

Develop on `nightly`; `main` is the release branch. Both branches run the regression
checks and build the documentation. Only `main` deploys the public documentation
and can publish packages. Merge completed, checked changes into `main` promptly.
A merge without a version increase does not publish a package.

Commit as Einar Olafsson using `einar.olafsson@gmail.com`, without additional
co-author trailers. Use a clear subject and a body describing the reason for the
change, resulting behaviour, and relevant validation. Automated publishing creates
no commits and does not add a code contributor.

## Local checks

Use a virtual environment and install the development and documentation tools:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,docs]"
python scripts/release.py check
QT_QPA_PLATFORM=offscreen pytest tests/test_release.py tests/test_wheel.py tests/test_packaging.py tests/test_readme.py tests/test_docstrings.py tests/test_paths.py tests/test_app_smoke.py tests/test_app_polish.py tests/test_discover_cli.py tests/test_importer.py -q
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

There is one PyPI distribution: **`starplast`**. It contains the Python modules,
console commands, built gene data, and icons. `starplast[gpu]` enables optional
CUDA dependencies and `starplast[ingest]` enables coverage-file imports; both are
extras of the same project.

The wheel includes cached data for offline browsing and excludes local saved
embeddings. Original source datasets are not packaged; see the [dataset catalogue](datasets.md).

## One-time PyPI setup

Sign in to **Einar's own PyPI account** at
[PyPI publishing](https://pypi.org/manage/account/publishing/) and create
a pending Trusted Publisher for **`starplast`** with these values:

| Field | Value |
|---|---|
| PyPI project name | `starplast` |
| Owner | `EinarOlafsson` |
| Repository | `starplast` |
| Workflow filename | `release.yml` |
| Environment | `pypi` |

For an existing project, add the same publisher in that project's Publishing
settings. See [PyPI's setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).
The workflow exchanges its GitHub identity for short-lived upload credentials;
no permanent PyPI token is required in GitHub secrets.
The account that registers a pending publisher owns the resulting PyPI project.
The `Owner` field above identifies the GitHub repository owner; it is not a PyPI
username. Do not create the project under a different account.

Create the GitHub environment `pypi` and restrict deployments to the `main` branch.
To publish automatically, leave required reviewers unset. Under repository
**Settings → Pages**, choose **GitHub Actions** as the source. The documentation
workflow builds pull requests and deploys pushes to `main`.

## Publish a version

```bash
git switch nightly
git pull --ff-only origin nightly
python scripts/release.py bump 0.43.0
python scripts/release.py check
```

Use a version greater than the current one. The bump command updates `pyproject.toml` and `starplast.__version__`.
Update `CHANGELOG.md`, run the release checks, commit with a descriptive message,
and push to `nightly`. Merge `nightly` into `main` through a pull request or locally:

```bash
git push origin nightly
git switch main
git pull --ff-only origin main
git merge --ff-only nightly
git push origin main
git switch nightly
```

If the branches have diverged, resolve the merge on `nightly` or use a pull request;
do not force-push `main`. A separate version-bump commit on `main` also triggers
the release, but preparing it on `nightly` keeps the normal development flow intact.

The release workflow compares the version with the repository state before the
push. An unchanged version does not publish; a downgrade fails. A version increase
runs the checks, builds a wheel and source distribution, verifies their contents,
and publishes `starplast`. After a successful PyPI upload, it creates
a GitHub release and `v<version>` tag at the exact commit that was built, with
generated release notes and the distribution files attached. Prerelease versions
are marked as prereleases on GitHub. Existing releases and their assets are left
intact on retries. Tags do not need to be created manually.

For the first upload or to retry a partially completed upload, manually run
**Publish Python packages** on `main`. Existing files are skipped; PyPI does not
allow an uploaded version to be replaced. Use a new version for changed code.
If publication failed before any upload, rerunning the original job preserves the
same tested artifacts.

## Inspect release artifacts

```bash
python -m build
python scripts/check_wheel.py dist
python -m twine check --strict dist/*
```

Start with an empty `dist/` directory so files from older versions are not included.
The wheel check requires both organism caches, the application icon and diagrams,
and the compressed sequence tables. It rejects saved embeddings, unconditional
CUDA dependencies, and files above PyPI's default upload limit.

The desktop installers in `packaging/` are separate from PyPI. Their existing
workflow runs manually or on a `v*` tag.
