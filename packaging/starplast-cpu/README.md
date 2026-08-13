# starplast-cpu

`pip install starplast-cpu` installs [starplast](https://github.com/EinarOlafsson/starplast) and
nothing else.

The program has never required a GPU. This name exists so that the pair reads symmetrically beside
`starplast-gpu`, and so that "not the two gigabytes of CUDA wheels, thank you" is something you can
say in the install command.

To add GPU support later, on the machine that will run it:

    starplast-install-gpu        # detects the driver, picks the CUDA set, shows the command first
