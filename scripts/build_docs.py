"""Build the user guides and pdoc API reference into docs/site."""
from __future__ import annotations

import html
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PAGES = {"index": "docs/index.md", "guide": "docs/guide.md", "API": "docs/API.md",
         "datasets": "docs/datasets.md", "releases": "docs/releases.md",
         "repository-review": "docs/repository-review.md",
         "scientific-roadmap": "docs/scientific-roadmap.md", "changelog": "CHANGELOG.md"}


def main():
    """Generate import-based API pages and render the Markdown guides beside them."""
    import markdown2
    output = ROOT / "docs" / "site"
    output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", PYQTGRAPH_QT_LIB="PyQt6")
    subprocess.run([sys.executable, "-m", "pdoc", "--output-directory", str(output / "api"),
                    "--docformat", "markdown", "--no-show-source",
                    "--edit-url", "starplast=https://github.com/EinarOlafsson/starplast/blob/main/starplast/",
                    "starplast"],
                   cwd=ROOT, env=env, check=True)
    shutil.copytree(ROOT / "docs" / "assets", output / "assets", dirs_exist_ok=True)
    shutil.copytree(ROOT / "docs" / "screenshots", output / "screenshots", dirs_exist_ok=True)
    shutil.copy2(ROOT / "starplast/data/icons/starplast.svg", output / "assets/icon.svg")
    css = '''body{margin:0;background:#f6f8fa;color:#192c3e;font:17px/1.65 system-ui,sans-serif}
    main{max-width:980px;margin:auto;padding:32px}nav{display:flex;gap:20px;flex-wrap:wrap}
    nav a{font-weight:600}a{color:#066b72}h1,h2,h3{line-height:1.2;color:#102a43}
    h1{font-size:2.5rem;margin-top:42px}h2{margin-top:40px}img{max-width:100%;height:auto}
    pre{padding:20px;background:#102a43;color:#e7f5f5;overflow-x:auto;border-radius:8px}
    code{font-size:.9em}table{display:block;overflow-x:auto;border-collapse:collapse;width:100%}
    td,th{padding:12px;border-bottom:1px solid #d7e1e6;text-align:left;vertical-align:top}
    th{background:#e5eff1}footer{margin-top:48px;color:#526879}'''
    (output / "style.css").write_text(css, encoding="utf-8")
    for page, source in PAGES.items():
        text = (ROOT / source).read_text(encoding="utf-8")
        # Relative guide links resolve both in the repository and in the built site.
        text = re.sub(r'\]\((?:docs/)?([\w-]+)\.md(#[^)]*)?\)',
                      lambda m: ("](" + ((m[1] if m[1] != "CHANGELOG" else "changelog") + ".html"
                                 if m[1] in {*PAGES, "CHANGELOG"}
                                 else "https://github.com/EinarOlafsson/starplast/blob/main/" + m[1] + ".md")
                                 + (m[2] or "") + ")"),
                      text)
        text = text.replace("](docs/assets/", "](assets/")
        body = markdown2.markdown(text, extras=["fenced-code-blocks", "tables", "header-ids"])
        title = html.escape(text.splitlines()[0].lstrip("# "))
        rendered = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} · Starplast</title><link rel="icon" href="assets/icon.svg" type="image/svg+xml">
<link rel="stylesheet" href="style.css"></head><body><main>
<nav aria-label="Main"><a href="index.html">Starplast</a><a href="guide.html">Guide</a>
<a href="API.html">API</a><a href="datasets.html">Datasets</a><a href="changelog.html">Changes</a>
<a href="https://github.com/EinarOlafsson/starplast">GitHub</a></nav>
{body}<footer>Starplast · Einar Olafsson</footer></main></body></html>'''
        (output / f"{page}.html").write_text(rendered, encoding="utf-8")
    print(f"Documentation: {output / 'index.html'}")


if __name__ == "__main__":
    main()
