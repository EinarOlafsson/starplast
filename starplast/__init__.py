"""Explore gene evidence in Toxoplasma gondii and Plasmodium falciparum.

Use the analysis modules from Python, ``starplast-discover`` for batch searches,
or ``starplast`` for the desktop browser. The user guide and API examples are
available at https://einarolafsson.github.io/starplast/.

The main modules load on first use, so ``import starplast`` stays light and never imports Qt::

    import starplast
    starplast.strategies.overview("Tg")          # every strategy, its task, grade and skill
    starplast.strategies.test("stacking", target="compartment").card()
    starplast.scorecard.glossary("ranking")      # what each metric means
"""
__version__ = "0.46.0"

#: Submodules reachable as attributes of the package, imported when first touched.
_LAZY = ("strategies", "scorecard", "techniques", "calibration", "graphspace", "deposits",
         "datasets", "slots", "paths", "search", "leakage")


def __getattr__(name: str):
    if name in _LAZY:
        import importlib
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module 'starplast' has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(_LAZY))
