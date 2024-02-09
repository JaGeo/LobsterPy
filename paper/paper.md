---
title: 'LobsterPy: A package to automatically analyze LOBSTER runs'
tags:
  - Python
  - Automation
  - Bonding analysis
  - Machine learning
authors:
  - name: Aakash Ashok Naik
    orcid: 0000-0002-6071-6786
    affiliation: "1 , 2"
  - name: Katharina Ueltzen
    orcid: 0009-0003-2967-1182
    affiliation: 1
  - name: Christina Ertural
    orcid: 0000-0002-7696-5824
    affiliation: 1
  - name: Adam J. Jackson
    orcid: 0000-0001-5272-6530
    affiliation: 3
  - name: Janine George
    orcid: 0000-0001-8907-0336
    affiliation: "1, 2"
affiliations:
 - name: Federal Institute for Materials Research and Testing, Department Materials Chemistry, Berlin, 12205, Germany
   index: 1
 - name: Friedrich Schiller University Jena, Institute of Condensed Matter Theory and Solid-State Optics, Jena, 07743, Germany
   index: 2
 - name: Scientific Computing Department, Science and Technology Facilities Council, Rutherford Appleton Laboratory, Didcot, 0X11 0QX, UK
   index: 3
date: August 2023
bibliography: paper.bib

---
# Summary
The LOBSTER[@deringer2011crystal;@maintz2013analytic;@maintz2016lobster;@nelson2020lobster] software aids in extracting quantum-chemical bonding information from materials by projecting the
plane-wave based wave functions from density functional theory (DFT) onto an atomic orbital basis. [LobsterEnv](https://github.com/materialsproject/pymatgen/blob/master/pymatgen/io/lobster/lobsterenv.py),
a module implemented in pymatgen[@ong2013python] by some of the authors of this package, facilitates the use of quantum-chemical bonding
information obtained from LOBSTER calculations to identify neighbors and coordination environments. _LobsterPy_ is a Python package that offers a set of convenient tools
to further analyze and summarize the LobsterEnv outputs in the form of JSONs that are easy to interpret and process. These tools enable the
estimation of (anti) bonding contributions, generation of textual descriptions, and visualization of LOBSTER computation results.  Since its first release, both _LobsterPy_ and _LobsterEnv_ capabilities
have been extended significantly. Unlike earlier versions, which could only automatically analyze Crystal Orbital Hamilton Populations (COHPs)[@dronskowski1993crystal],
both can now also analyze Crystal Orbital Overlap Populations (COOP)[@hughbanks1983chains] and Crystal Orbital Bond Index (COBI)[@mueller2021crystal].
Extracting the information about the most important orbitals contributing to the bonds is optional, and users can enable it as needed.
Additionally, bonding-based features for machine-learning (ML) studies can be engineered via the sub-packages "featurize" and "structuregraphs".
Alongside its Python interface, it also provides an easy-to-use command line interface (CLI) that runs automatic analysis of the
computations and generates a summary of results and publication-ready figures.

_LobsterPy_ has been used to produce the results in [@ngo2023dft; @chen2024insights; @naik2023quantumchemical] and is also part of
[@atomate2] bonding analysis workflow for generating bonding analysis data in a format compatible with the Materials Project[@materialsproject] API.

# Statement of need
Although the notion of "bonds" might seem unusual from a physicist's point of view, chemists have been employing it routinely to
explain various chemical phenomena and materials properties.[@hoffmann1987chemistry; @burdett1995chemical; @das2023strong; @ertural2022first; @hu2023mechanism; @dronskowski2023chemical] With the recent advances in
automation frameworks for high-throughput computational investigations, bonding analysis for thousands of crystalline materials
can be performed with few lines of code.[@george2022automated] This automation helps reduce the common mistakes inexperienced
users make while performing bonding analysis. However, it is also essential to systematically generate inputs and post-process
the output files consistently to have reliable and reproducible results. Furthermore, transforming the data from these high-throughput
bonding analysis calculations into a format suitable for ML studies should benefit data-driven material science research.
_LobsterPy_ aims to fulfill this need.

# Features
- Generate summarized bonding analysis JSONs and text descriptions based on COHPs (ICOHPs), COBIs (ICOBIs), and COOPs (ICOOPs)
- Generate static and interactive plots of the most relevant COHPs, COBIs, and COOPs
- Customizable plotters for visualization of COHPs (ICOHPs), COBIs (ICOBIs), COOPs (ICOOPs) and DOS
- Benchmark LOBSTER calculation quality and generate corresponding JSONs and text descriptions
- Create inputs for LOBSTER calculations from VASP files
- Extract features from LOBSTER calculation files to be used for ML studies
- Perform automatic bonding analysis and plotting via inherent command line interface app.

# Similar and Related Software
LobsterPy can be seen to be similar in spirit to sumo [@Ganose2018], as both provide Python tools to analyze and visualize data related to the electronic structure that are based on ab initio calculations. Other software packages that enable visualizing results specifically from the LOBSTER software are wxDragon[@wxdragon] and [@abipy].LobsterPy differs from these two packages by providing further analysis of the calculations, interpretable text summaries, and featurizers for ML studies besides plotting the data.

# Availability
_LobsterPy_ can also be found on [PyPI](https://pypi.org/project/lobsterpy/). Detailed software documentation,
including [installation instructions](https://jageo.github.io/LobsterPy/installation/index.html) and
[implementation details](https://jageo.github.io/LobsterPy/fundamentals/index.html) are provided. The package
also includes [tutorials](https://jageo.github.io/LobsterPy/tutorial/index.html) illustrating all the basic and
advanced functionalities.

# Acknowledgements
The authors would like to acknowledge the Gauss Centre for Super
computing e.V. (www.gauss-centre.eu) for funding this project by
providing generous computing time on the GCS Supercomputer
SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de)
(project pn73da) that enabled rigorous testing of this
package on a diverse set of compounds. The authors thank Jonas Grandel for reviewing the docstrings and testing package functionalities
and tutorials. The authors would also like to acknowledge the maintainers of pymatgen and LOBSTER program developers.

# References
