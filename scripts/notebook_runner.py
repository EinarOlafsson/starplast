"""Write an executed Jupyter notebook from Python, without Jupyter.

The project rule is that every data download, analysis and figure is kept as an annotated notebook
in the tree. This builds one by running each code cell in a shared namespace, capturing what it
prints and the value of its last expression, and writing the .ipynb with those outputs -- so the
notebook on disk is what actually ran, not a template someone forgot to execute.

    nb = ExecutedNotebook("title")
    nb.md("## A heading", "Why this step exists.")
    nb.code("x = 1 + 1", "x")          # lines joined; the last expression becomes the output
    nb.write("notebooks/example.ipynb")
"""
from __future__ import annotations

import ast
import base64
import contextlib
import io
import json
import os
import time


class ExecutedNotebook:
    def __init__(self, title: str, namespace: dict | None = None):
        self.title = title
        self.ns = {} if namespace is None else namespace
        self.cells = []
        self.count = 0

    def md(self, *paragraphs: str):
        text = "\n\n".join(paragraphs)
        self.cells.append({"cell_type": "markdown", "metadata": {},
                           "source": text.splitlines(keepends=True)})

    def code(self, *lines: str):
        import pandas as pd
        source = "\n".join(lines)
        self.count += 1
        tree = ast.parse(source)
        body, last = source, None
        if tree.body and isinstance(tree.body[-1], ast.Expr):
            start = tree.body[-1].lineno - 1
            src_lines = source.splitlines()
            body, last = "\n".join(src_lines[:start]), "\n".join(src_lines[start:])
        buf = io.StringIO()
        t0 = time.monotonic()
        value = None
        with contextlib.redirect_stdout(buf):
            exec(compile(body, f"<cell {self.count}>", "exec"), self.ns)
            if last is not None:
                value = eval(compile(last, f"<cell {self.count}>", "eval"), self.ns)
        outputs = []
        if buf.getvalue():
            outputs.append({"output_type": "stream", "name": "stdout", "text": buf.getvalue()})
        if value is not None:
            data = {}
            if hasattr(value, "savefig"):
                png = io.BytesIO()
                value.savefig(png, format="png", dpi=90, bbox_inches="tight")
                data = {"image/png": base64.b64encode(png.getvalue()).decode(),
                        "text/plain": "<figure>"}
                import matplotlib.pyplot as plt
                plt.close(value)
                outputs.append({"output_type": "display_data", "metadata": {}, "data": data})
            else:
                if isinstance(value, (pd.DataFrame, pd.Series)):
                    shown = value.head(20) if len(value) > 20 else value
                    data = {"text/plain": shown.to_string()}
                    if isinstance(shown, pd.DataFrame):
                        data["text/html"] = shown.to_html(float_format=lambda x: f"{x:.4g}")
                else:
                    data = {"text/plain": repr(value)}
                outputs.append({"output_type": "execute_result", "execution_count": self.count,
                                "metadata": {}, "data": data})
        self.cells.append({"cell_type": "code", "execution_count": self.count,
                           "metadata": {"seconds": round(time.monotonic() - t0, 2)},
                           "source": source.splitlines(keepends=True), "outputs": outputs})
        return value

    def write(self, path: str):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        nb = {"cells": [{"cell_type": "markdown", "metadata": {},
                         "source": [f"# {self.title}\n"]}] + self.cells,
              "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                          "name": "python3"},
                           "language_info": {"name": "python"}},
              "nbformat": 4, "nbformat_minor": 5}
        with open(path, "w") as fh:
            json.dump(nb, fh, indent=1)
        return path
