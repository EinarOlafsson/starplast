# Space packs

Large organism spaces can live outside the application wheel. The pack framework is implemented;
human/mouse gene-space builders, published pack URLs and the download menu remain pending. The
human/mouse protein reference tables already shipped with Starplast are separate from these future
gene spaces.

A pack is a ZIP containing flat data files and `manifest.json`. The manifest records the organism,
version, node/graph filenames, each file's SHA256 and byte size, and a license code for every node
column. Only reviewed open (`O`) or share-alike (`SA`) columns can be packaged for distribution.
The builder's license declarations must come from the source records; hashes do not establish
redistribution rights or publisher identity.

Validate and install a local pack:

```bash
python -m starplast.packs verify /path/to/pack.zip --sha256 TRUSTED_SHA256
python -m starplast.packs install /path/to/pack.zip --sha256 TRUSTED_SHA256 --organism CODE
python -m starplast.packs installed CODE
```

For a published HTTPS URL, use `download URL --sha256 TRUSTED_SHA256 --organism CODE`.
Download requires a checksum obtained independently from the archive. Local installation can omit
the archive checksum when the local input is trusted; member checksums are always verified.

Installed versions live at `user_cache_dir()/spaces/CODE/VERSION/`. `STARPLAST_STATE` overrides the
user cache root. Installation checks the entire archive before activation, refuses unsafe members,
and replaces `active.json` atomically. Previous versions remain available. Reinstalling identical
data is harmless; different data under the same version is refused. Reinstall an older archive to
select it again. No imports or path lookups download anything.

Prepare a pack for a registered organism:

```bash
python scripts/build_space.py --organism CODE --data /path/to/prepared/tables \
  --version VERSION --licenses /path/to/column-licenses.json --output /path/to/pack.zip
```

The license file is a JSON mapping of every node column to its license code. The shared builder
requires unique canonical gene IDs, matching graph order, finite coordinates and valid edge
indices. It compares with the installed table by default, or with `--previous TABLE.parquet`, and
refuses lost identifiers, columns or previously measured cells. This packages prepared data; it
does not acquire or derive new measurements.

Python callers use `starplast.spaces.build_pack`, `starplast.packs.install` and
`starplast.organisms.nodes_path` / `graph_path`. A registry entry with `distribution="pack"` resolves
only its active installed pack. Wheel spaces retain the existing `STARPLAST_CACHE` resolution.
