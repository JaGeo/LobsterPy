![CI Status](https://github.com/JaGeo/LobsterPy/actions/workflows/python-package.yml/badge.svg) [![Documentation Status](https://readthedocs.org/projects/lobsterpy/badge/?version=latest)](https://lobsterpy.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/lobsterpy.svg)](https://badge.fury.io/py/lobsterpy) [![PyPI downloads](https://img.shields.io/pypi/dm/lobsterpy?style=flat&color=blue&label=pypi%20downloads)](https://pypi.org/project/lobsterpy) [![Downloads](https://pepy.tech/badge/lobsterpy)](https://pepy.tech/project/lobsterpy)   [![DOI](https://zenodo.org/badge/343384088.svg)](https://zenodo.org/badge/latestdoi/343384088)

# Getting started
<img src="https://raw.githubusercontent.com/JaGeo/LobsterPy/main/LobsterPyLogo.png" alt="LobsterPy Logo which consists of a green Python and a red Lobster" width="200"/>

This is a package that enables automatic plotting and summaries of Lobster outputs. You can download Lobster on [http://www.cohp.de](http://www.cohp.de). Currently, only VASP/Lobster computations are supported.

Please note that LobsterPy relies on your Lobster outputs. Thus, please make sure that the outputs have enough information for our (automatic) analysis.

![LobsterPyAnimation](https://github.com/JaGeo/LobsterPy/assets/22094846/8f06b84c-db6d-414c-8590-aa04c957c728)


## Installation

You can now use ``pip install lobsterpy`` to install it.

You can also pip install the package in development mode by writing ``pip install -e .``. It will then use setup.py to install the package. One requirement of this package is [pymatgen](https://github.com/materialsproject/pymatgen).

## Basic usage

* **Automatic analysis and plotting of COHPs:**

    <img src="https://github.com/JaGeo/LobsterPy/assets/22094846/6587e752-6ea4-4358-a763-3633d5a21869" alt="Output Automatic Analysis" width="300"/>

    
    You can use ``lobsterpy description`` for an automated analysis of COHPs for relevant cation-anion bonds or ``lobsterpy automatic-plot`` to plot the results automatically. It will evaluate all COHPs with ICOHP values down to 10% of the strongest ICOHP. You can enforce an analysis of all bonds by using ``lobsterpy automatic-plot --allbonds`` . Currently, the computed Mulliken charges will be used to determine cations and anions. If no ``CHARGE.lobster`` is available, the algorithm will fall back to the BondValence analysis from pymatgen. Please be aware that LobsterPy can only analyze bonds that have been included in the initial Lobster computation. Thus, please use the cohpgenerator within Lobster.

    It is also possible to start this automatic analysis from Python script. See "examples" for scripts.


* **Command line plotter**:

    We included options to plot COHPs/COBIs/COOPs from the command line.
    ``lobsterpy plot 1 2`` will plot COHPs of the first and second bond from ``COHPCAR.lobster``. It is possible to sum or integrate the COHPs as well (``--summed``, ``--integrated``). You can switch to COBIs or COOPs by using ``--cobis`` or ``--coops``, respectively.

* **Other command line tools**:

    ``lobsterpy create-inputs`` will create standard inputs based on existing POSCAR, POTCAR, INCAR files. It will allow to test for different basis sets that are available in Lobster. Currently only available for PBE_54 POTCARs.



* **Further help?**

    You can get further information by using ``lobsterpy --help`` and also by typing ``lobsterpy description --help``, ``lobsterpy automatic-plot --help``, ``lobsterpy plot --help``

## How to cite?
Please cite our paper: J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier, **ChemPlusChem**, [https://doi.org/10.1002/cplu.202200123](https://doi.org/10.1002/cplu.202200123).
Please cite [pymatgen](https://github.com/materialsproject/pymatgen), [Lobster](https://www.cohp.de), and [ChemEnv](https://doi.org/10.1107/S2052520620007994) correctly as well.

## LobsterPy as part of an atomate2 workflow
![LobsterWorkflow](https://github.com/JaGeo/LobsterPy/assets/22094846/337615ac-542e-446c-bc63-fb5946b16544)

We have now also included the automatic analysis into a fully automatic workflow using VASP and Lobster in [atomate2](https://github.com/materialsproject/atomate2). More documentation and information will follow soon.

## Future plans:
* Include automatic plotting for COBIs/COOPs
* Include orbitals into automatic plotting

## Contributions
* Contributions and suggestions for features are also welcome. Please write an Issue to describe your potential contribution or feature request.
* We are planning to submit a paper for the code LobsterPy when more features have been added (~ mid of 2023). Major contributors will of course have the chance to be co-authors. Please talk to us if you are interested in contributing :).

## Acknowledgements
The development of the program has been supported by a computing time grant. We gratefully acknowledge the Gauss Centre for Supercomputing e.V.(www.gauss-centre.eu) for funding this project by providing computing time on the GCS Supercomputer SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de) (project pn73da).
