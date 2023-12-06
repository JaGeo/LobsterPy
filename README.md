![CI Status](https://github.com/JaGeo/LobsterPy/actions/workflows/python-package.yml/badge.svg) [![build-docs](https://github.com/JaGeo/LobsterPy/actions/workflows/docs.yml/badge.svg)](https://jageo.github.io/LobsterPy/) [![PyPI version](https://badge.fury.io/py/lobsterpy.svg)](https://badge.fury.io/py/lobsterpy) [![PyPI downloads](https://img.shields.io/pypi/dm/lobsterpy?style=flat&color=blue&label=pypi%20downloads)](https://pypi.org/project/lobsterpy) [![Downloads](https://pepy.tech/badge/lobsterpy)](https://pepy.tech/project/lobsterpy)   [![DOI](https://zenodo.org/badge/343384088.svg)](https://zenodo.org/badge/latestdoi/343384088)

# Getting started
<img src="https://raw.githubusercontent.com/JaGeo/LobsterPy/main/LobsterPyLogo.png" alt="LobsterPy Logo which consists of a green Python and a red Lobster" width="200"/>

LobsterPy is a package that enables automatic analysis of LOBTSER outputs to get summarized bonding information and relevant bond plots. Additionally, one can also generate features for machine learning studies from LOBSTER outputs. One can download LOBSTER from [http://www.cohp.de](http://www.cohp.de).

<div style="border: 1px solid #52a5ab; padding: 5px; position: relative;">
    <div style="background-color: #52a5ab; color: #ffffff; padding: 0px; position: absolute; top: 0; left: 0; right: 0; text-align: center;">
        <strong>Important</strong>
    </div>
<br>

Recently released [LOBSTER 5.0](https://schmeling.ac.rwth-aachen.de/cohp/index.php?menuID=6) now generates `POSCAR.lobster` for any kind of LOBSTER calculation by default (This file has same format as the POSCAR from VASP). Thus, LobsterPy in principle, now supports usage with **all** DFT codes supported by LOBSTER and is **no** longer limited to `VASP`. Almost all of the core functionalities of LobsterPy could be used. The user must use `POSCAR.lobster` for `path_to_poscar` and `-fstruct` argument in python and cli interface, respectively.

The only functionality limited to VASP is DOS comparisons and basis set analysis in the `calc_quality_summary` method of the `Analysis` class, as it relies on VASP output files, namely `vasprun.xml` and `POTCAR`.
</div>

Please note that LobsterPy relies on the LOBSTER computation output files. Thus, it will be only able to analyze data that has been computed in the LOBSTER run.

![LobsterPyAnimation](https://github.com/JaGeo/LobsterPy/assets/22094846/8f06b84c-db6d-414c-8590-aa04c957c728)


## Installation
### Standard installation
Install using ``pip install lobsterpy``

### Installation with featurizer
Install using ``pip install lobsterpy[featurizer]``


## Basic usage

* **Automatic analysis and plotting of COHPs / COBIS / COOPs:**

    <img src="https://github.com/JaGeo/LobsterPy/assets/22094846/6587e752-6ea4-4358-a763-3633d5a21869" alt="Output Automatic Analysis" width="300"/>

You can use ``lobsterpy description`` for an automated analysis of COHPs for relevant cation-anion bonds or ``lobsterpy automatic-plot`` to plot the results automatically.
It will evaluate all COHPs with ICOHP values down to 10% of the strongest ICOHP.
You can enforce an analysis of all bonds by using ``lobsterpy automatic-plot --allbonds``.
You can also switch the automatic analysis to use the ICOBIs or ICOOPs. You need to add `--cobis` or `--coops` along with the mentioned commands
for e.g.like  ``lobsterpy description --cobis``

An interactive plotter is available via ``lobsterpy automatic-plot-ia``.

Currently, the computed Mulliken charges will be used to determine cations and anions. If no ``CHARGE.lobster`` is available, the algorithm will fall back to the BondValence analysis from pymatgen.

*Please be aware that LobsterPy can only analyze bonds that have been included in the initial Lobster computation. Thus, please use the cohpgenerator within Lobster.*


It is also possible to start this automatic analysis from a Python script. See "examples" for scripts.

* **Plotting DOS from LOBSTER computations:**

  To plot densities of states obtained from LOBSTER use ``lobsterpy plot-dos``.


* **Generic COHP/ COOP / COBI plotter**:

  We included options to plot COHPs/COBIs/COOPs from the command line.
``lobsterpy plot 1 2`` will plot COHPs of the first and second bond from ``COHPCAR.lobster``. It is possible to sum or integrate the COHPs as well (``--summed``, ``--integrated``). You can switch to COBIs or COOPs by using ``--cobis`` or ``--coops``, respectively.

* **Other command line tools**:

    ``lobsterpy create-inputs`` will create standard inputs based on existing POSCAR, POTCAR, INCAR files. It will allow to test for different basis sets that are available in Lobster. Currently only available for PBE_54 POTCARs.


* **Further help?**

    You can get further information by using ``lobsterpy --help`` and also by typing ``lobsterpy description --help``,
``lobsterpy automatic-plot --help``, ``lobsterpy plot --help``.

## Comprehensive documentation
* Checkout the [documentation and tutorials](https://jageo.github.io/LobsterPy/) for more details.

## How to cite?
Please cite our paper: J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier, **ChemPlusChem**, [https://doi.org/10.1002/cplu.202200123](https://doi.org/10.1002/cplu.202200123).
Please cite [pymatgen](https://github.com/materialsproject/pymatgen), [Lobster](https://www.cohp.de), and [ChemEnv](https://doi.org/10.1107/S2052520620007994) correctly as well.

## LobsterPy is now a part of an atomate2 workflow
![LobsterWorkflow](https://github.com/JaGeo/LobsterPy/assets/22094846/337615ac-542e-446c-bc63-fb5946b16544)

We have now also included the automatic analysis into a fully automatic workflow using VASP and Lobster in [atomate2](https://github.com/materialsproject/atomate2). More documentation and information will follow soon.


## Acknowledgements
The development of the program has been supported by a computing time grant. We gratefully acknowledge the Gauss Centre for Supercomputing e.V.(www.gauss-centre.eu) for funding this project by providing computing time on the GCS Supercomputer SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de) (project pn73da).
