# Starplast

Explore gene evidence and screen results in *Toxoplasma gondii* and
*Plasmodium falciparum*.

```bash
pip install starplast
starplast
```

Requires Python 3.10+, a display, and OpenGL. Built gene tables and graphs are
included for offline browsing. Analysis runs on the CPU by default.
Install `starplast[gpu]` for optional CUDA 12 acceleration on Linux x86_64.

[Documentation](https://einarolafsson.github.io/starplast/) ·
[Source code](https://github.com/EinarOlafsson/starplast)

This package installs the matching `starplast-core` distribution, which contains
the code and data. The `gpu` and `ingest` extras forward to that distribution.
