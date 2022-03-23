![CI Status](https://github.com/JaGeo/LobsterPy/actions/workflows/python-package.yml/badge.svg)  [![PyPI version](https://badge.fury.io/py/lobsterpy.svg)](https://badge.fury.io/py/lobsterpy) [![PyPI downloads](https://img.shields.io/pypi/dm/lobsterpy?style=flat&color=blue&label=pypi%20downloads)](https://pypi.org/project/lobsterpy) [![DOI](https://zenodo.org/badge/343384088.svg)](https://zenodo.org/badge/latestdoi/343384088)

# LobsterPy
<img src="LobsterPyLogo.png" alt="LobsterPy Logo" width="200"/>
This is a package that enables automatic plotting of Lobster outputs. You can download Lobster on [http://www.cohp.de](http://www.cohp.de). Currently, only VASP/Lobster computations are supported.

## Installation

You can now use ``pip install lobsterpy`` to install it.

You can also pip install the package in development mode by writing ``pip install -e .``. It will then use setup.py to install the package. One requirement of this package is [pymatgen](https://github.com/materialsproject/pymatgen).

## Basic usage

* **Automatic analysis and plotting of COHPs:**
    
    You can use ``lobsterpy description`` for an automated analysis of COHPs for relevant cation-anion bonds or ``lobsterpy automatic-plot`` to plot the results automatically. It will evaluate all COHPs with ICOHP values down to 10% of the strongest ICOHP. You can enforce an analysis of all bonds by using ``lobsterpy automatic-plot --allbonds`` . Currently, the computed Mulliken charges will be used to determine cations and anions. If no ``CHARGE.lobster`` is available, the algorithm will fall back to the BondValence analysis from pymatgen.
  
    It is also possible to start this automatic analysis from Python script. See "examples" for scripts.

  
* **Command line plotter**:
    
    We included options to plot COHPs/COBIs/COOPs from the command line.
    ``lobsterpy plot 1 2`` will plot COHPs of the first and second bond from ``COHPCAR.lobster``. It is possible to sum or integrate the COHPs as well (``--summed``, ``--integrated``). You can switch to COBIs or COOPs by using ``--cobis`` or ``--coops``, respectively.


* **Further help?**
  
    You can get further information by using ``lobsterpy --help`` and also by typing ``lobsterpy description --help``, ``lobsterpy automatic-plot --help``, ``lobsterpy plot --help``


## License
Lobsterpy is released under a BSD 3-Clause "New" or "Revised" License. 


## How to cite?
Please cite our preprint: [https://doi.org/10.26434/chemrxiv-2022-2v424](https://doi.org/10.26434/chemrxiv-2022-2v424). 
Please cite [pymatgen](https://github.com/materialsproject/pymatgen), [Lobster](https://www.cohp.de), and [ChemEnv](https://doi.org/10.1107/S2052520620007994) correctly as well.


## Future plans:
* Include automatic plotting for COBIs/COOPs
* Include orbitals into automatic plotting
* Include more documentation
* Include lobsterin generation

## Contributions
* Contributions and suggestions for features are also welcome. Please write an Issue to describe your potential contribution or feature request.
