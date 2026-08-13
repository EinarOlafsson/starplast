# starplast

`pip install starplast` installs the program **and** the GPU stack it can use, where this platform
has wheels for it: cuML for UMAP and HDBSCAN, CuPy for the array work. On macOS and Windows, where
RAPIDS publishes nothing, it installs the program alone rather than failing.

    pip install starplast          the program, GPU where possible
    pip install starplast-cpu      the program alone, no CUDA wheels
    starplast-install-gpu          add GPU support later, choosing the CUDA set from the driver

The code itself is the `starplast-core` distribution; the import name is `starplast` either way.
