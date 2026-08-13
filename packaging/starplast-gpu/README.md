# starplast-gpu

`pip install starplast-gpu` installs [starplast](https://github.com/EinarOlafsson/starplast) together
with the CUDA 12 stack it can use: **cuML** for UMAP and HDBSCAN, **CuPy** for the array work.

It contains no code of its own. Everything it installs is `starplast[gpu]`.

Nothing about starplast requires it. With no GPU stack present the program runs exactly as it does
now — `Preferences ▸ compute` reports what it found, and the switch has nothing to turn on.
