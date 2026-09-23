"""Refresh the dataset catalogue from the application registry."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    """Replace the generated catalogue while retaining its introduction."""
    from starplast.datasets import readme_table
    target = ROOT / "docs" / "datasets.md"
    intro = target.read_text(encoding="utf-8").split("### ", 1)[0]
    target.write_text(intro + readme_table().strip() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
