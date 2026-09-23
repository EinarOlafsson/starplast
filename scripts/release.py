"""Synchronize release versions and decide whether a push should publish packages."""
from __future__ import annotations

import argparse
import ast
from pathlib import Path
import re
import subprocess

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[1]
PROJECTS = (Path("pyproject.toml"),)


def check(root: Path = ROOT) -> str:
    """Return the version after checking package identity and the runtime declaration."""
    projects = [tomllib.loads((root / p).read_text(encoding="utf-8"))["project"] for p in PROJECTS]
    version = projects[0]["version"]
    if str(Version(version)) != version:
        raise ValueError(f"Use a canonical PEP 440 version: {version}")
    for p in projects:
        if p["version"] != version:
            raise ValueError(f"{p['name']} has version {p['version']}; expected {version}")
    tree = ast.parse((root / "starplast/__init__.py").read_text(encoding="utf-8"))
    runtime = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == "__version__" for t in n.targets))
    if runtime != version:
        raise ValueError(f"starplast.__version__ is {runtime}; expected {version}")
    project = projects[0]
    if project["name"] != "starplast":
        raise ValueError("The only PyPI project must be starplast")
    requirements = project.get("dependencies", []) + [
        r for values in project.get("optional-dependencies", {}).values() for r in values]
    if any(r.lower().startswith("starplast") for r in requirements):
        raise ValueError("Starplast must not depend on another Starplast distribution")
    return version


def bump(version: str, root: Path = ROOT) -> str:
    """Update package and runtime declarations to a newer PEP 440 version."""
    old = check(root)
    if str(Version(version)) != version or Version(version) <= Version(old):
        raise ValueError(f"Version must be canonical and greater than {old}")
    changes = {}
    for path in PROJECTS:
        text = (root / path).read_text(encoding="utf-8")
        text = re.sub(r'(?m)^version = "[^"]+"$', f'version = "{version}"', text)
        changes[path] = text
    path = Path("starplast/__init__.py")
    changes[path] = re.sub(r'(?m)^__version__ = "[^"]+"$', f'__version__ = "{version}"',
                           (root / path).read_text(encoding="utf-8"))
    for path, text in changes.items():
        (root / path).write_text(text, encoding="utf-8")
    return check(root)


def changed(previous: str, root: Path = ROOT) -> bool:
    """Publish only a version increase relative to the state before a push; reject downgrades."""
    current = Version(check(root))
    if not re.fullmatch(r"[0-9a-fA-F]{40,64}", previous):
        raise ValueError("The previous revision must be a full Git object ID")
    if set(previous) == {"0"}:
        return True
    result = subprocess.run(["git", "show", f"{previous}:pyproject.toml"], cwd=root,
                            capture_output=True, text=True, check=True)
    old = Version(tomllib.loads(result.stdout)["project"]["version"])
    if current < old:
        raise ValueError(f"Refusing version downgrade: {old} -> {current}")
    return current > old


def main():
    """Check versions, bump them, or write publication eligibility to GitHub job outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check")
    update = sub.add_parser("bump")
    update.add_argument("version")
    detect = sub.add_parser("detect")
    detect.add_argument("--previous", default="")
    detect.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "bump":
        print(bump(args.version))
    else:
        version = check()
        if args.command == "detect":
            publish = changed(args.previous) if args.previous else True
            with args.output.open("a", encoding="utf-8") as stream:
                stream.write(f"version={version}\npublish={str(publish).lower()}\n")
                stream.write(f"prerelease={str(Version(version).is_prerelease).lower()}\n")
        print(version)


if __name__ == "__main__":
    main()
