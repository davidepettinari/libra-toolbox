# libra-toolbox

`libra-toolbox` is a Python package developed by the LIBRA team for neutron detection and tritium analysis. This toolbox provides various modules and functionalities to facilitate scientific research and data analysis for tritium breeding experiments.

## Installation

`libra-toolbox` relies on Python. You can install it from the official website: [python.org](https://www.python.org/downloads/).

Install OpenMC:

If you want to run OpenMC functions, [install OpenMC](https://docs.openmc.org/en/stable/quickinstall.html) with conda:

```
conda install -c conda-forge openmc>=0.14.0
```


To install `libra-toolbox`, use pip:

```bash
pip install git+https://github.com/LIBRA-project/libra-toolbox
```

## Documentation
The documentation for `libra-toolbox` is built using Sphinx and is [available online](https://libra-toolbox.readthedocs.io/en/latest/). To build the documentation locally, you can use the provided Makefile or make.bat script.

### Building Documentation
To build the documentation, navigate to the `docs` directory and run:

```
make html
```

or on Windows:

```
make.bat html
```

The generated documentation will be available in the `docs/_build/html` directory.