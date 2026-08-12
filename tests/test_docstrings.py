"""Every public module, class and function carries a docstring.

There will be an API before long, and the generated reference IS the documentation for this project
-- there is no second prose document to fall back on. A public name with no docstring is a hole in
that reference, and holes are much cheaper to prevent than to find later. 28 were missing when this
was first run, including the main Window class and the dataset record.

The bar is a docstring that exists and says something. Whether it says WHY rather than what is not
mechanically checkable, but it is the convention throughout: the module docstrings here carry the
reasoning, which is why pdoc's output is worth reading.
"""
import ast
import os

import pytest

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "starplast")
MODULES = sorted(f for f in os.listdir(ROOT) if f.endswith(".py"))


def _tree(name):
    return ast.parse(open(os.path.join(ROOT, name)).read())


@pytest.mark.parametrize("name", MODULES)
def test_every_module_explains_itself(name):
    doc = ast.get_docstring(_tree(name))
    assert doc and doc.strip(), f"{name} has no module docstring"


@pytest.mark.parametrize("name", MODULES)
def test_every_public_function_and_class_has_a_docstring(name):
    missing = []
    for node in _tree(name).body:
        public = not node.__dict__.get("name", "_").startswith("_")
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and public:
            if not ast.get_docstring(node):
                missing.append(f"{name}:{node.lineno} {node.name}")
    assert not missing, f"no docstring: {missing}"


@pytest.mark.parametrize("name", MODULES)
def test_public_methods_of_public_classes_have_docstrings(name):
    """The API surface people actually call, which is mostly methods rather than free functions."""
    missing = []
    for node in _tree(name).body:
        if not isinstance(node, ast.ClassDef) or node.name.startswith("_"):
            continue
        for sub in node.body:
            if not isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if sub.name.startswith("_") or ast.get_docstring(sub):
                continue
            # Qt event handlers and property getters are named by their framework, not by us, and a
            # docstring on `paintEvent` restating "paints the event" is noise.
            if sub.name.endswith("Event") or any(
                    getattr(d, "id", getattr(d, "attr", "")) == "property" for d in sub.decorator_list):
                continue
            missing.append(f"{name}:{sub.lineno} {node.name}.{sub.name}")
    assert not missing, f"no docstring: {missing}"
