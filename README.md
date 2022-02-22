![CI Status](https://github.com/JaGeo/LobsterPy/actions/workflows/python-package.yml/badge.svg)
 
# LobsterPy

This is a package that enables automatic plotting of Lobster outputs. You can download Lobster on [https://www.cohp.de](https://www.cohp.de). Currently, only VASP/Lobster computations are supported.

## Installation


You can pip install the package by writing ``pip install -e .`` in the main package. It will then use setup.py to install the package. One requirement of this package is [pymatgen](https://github.com/materialsproject/pymatgen).

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
A citation for lobsterpy will follow at a later stage. For now, please cite [pymatgen](https://github.com/materialsproject/pymatgen), [Lobster](https://www.cohp.de), and [ChemEnv](https://doi.org/10.1107/S2052520620007994) correctly.


## Future plans:
* Include automatic plotting for COBIs/COOPs
* Include orbitals into automatic plotting
* Include more documentation
