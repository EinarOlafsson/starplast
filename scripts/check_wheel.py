"""Check built package metadata, cache files, icons, and size before uploading."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
from pathlib import Path
import zipfile


def check_wheel(path: Path) -> None:
    """Reject incomplete wheels, accidental saved analyses, and oversized uploads."""
    if path.stat().st_size > 100_000_000:
        raise ValueError(f"{path.name} exceeds the default PyPI 100 MB file limit")
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        metadata = next(n for n in names if n.endswith(".dist-info/METADATA"))
        info = BytesParser().parsebytes(archive.read(metadata))
        if info["Name"] != "starplast":
            raise ValueError(f"Unexpected distribution: {info['Name']}")
        required = {"starplast/app.py", "starplast/data/nodes.parquet", "starplast/data/graph.npz",
                    "starplast/data/toxodb_identity.tsv", "starplast/data/pf_nodes.parquet",
                    "starplast/data/pf_graph.npz", "starplast/data/icons/starplast.svg",
                    "starplast/data/icons/Apicomplexa_cells.svg", "starplast/data/icons/Animal_cells.svg",
                    "starplast/data/toxodb_cds.tsv.gz", "starplast/data/plasmodb_cds.tsv.gz",
                    "starplast/workflows.py", "starplast/prediction.py", "starplast/evidence.py",
                    "starplast/data/af3_features.parquet", "starplast/data/af3_manifest.parquet",
                    "starplast/data/af3_features.json", "starplast/data/esm_features.parquet",
                    "starplast/data/esm_manifest.parquet", "starplast/data/pf_mentions.parquet"}
        missing = required - names
        if missing:
            raise ValueError(f"Incomplete wheel: {sorted(missing)}")
        if any(n.startswith("starplast/data/embeddings/") for n in names):
            raise ValueError("The wheel contains local saved embeddings")
        for requirement in info.get_all("Requires-Dist", []):
            if requirement.lower().startswith("starplast"):
                raise ValueError(f"Unexpected Starplast dependency: {requirement}")
            if requirement.lower().startswith(("cuml", "cupy")) and 'extra == "gpu"' not in requirement:
                raise ValueError(f"CUDA must be optional: {requirement}")
    print(f"Checked {path.name}: {path.stat().st_size / 1_000_000:.1f} MB")


def main():
    """Validate each wheel in a distribution directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path("dist"))
    directory = parser.parse_args().directory
    wheels = list(directory.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError(f"Expected one wheel, found {len(wheels)}")
    expected_sdist = wheels[0].name.split("-py", 1)[0] + ".tar.gz"
    if {p.name for p in directory.iterdir()} != {wheels[0].name, expected_sdist}:
        raise ValueError("Upload directory must contain only the matching starplast wheel and sdist")
    for wheel in wheels:
        check_wheel(wheel)


if __name__ == "__main__":
    main()
