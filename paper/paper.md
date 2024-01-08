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
 - name: Science and Technology Facilities Council, Didcot, Oxfordshire, GB
   index: 3
date: August 2023
bibliography: paper.bib

---
# Summary
The LOBSTER software aids in extracting quantum-chemical bonding information from materials. It does this by projecting the
plane-wave based wave functions from density functional theory (DFT) onto an atomic orbital basis.  _LobsterPy_ is a Python package
that provides convenient tools to systematically analyze, describe, and visualize such LOBSTER computations results.
Since its first release, its capabilities have been extended significantly. Unlike earlier versions, which could only
automatically analyze Crystal Orbital Hamilton Populations (COHPs)[@dronskowski1993crystal], _LobsterPy_ can now also analyze
Crystal Orbital Overlap Populations (COOP)[@hughbanks1983chains] and Crystal Orbital Bond Index (COBI)[@mueller2021crystal] to
extract summarized bonding information. The latter includes information on coordination environments, bond strengths, most relevant bonds,
bonding, and anti-bonding contributions. Optionally, users can further extract the most important orbitals contributing to the
relevant bonds. Additionally, bonding-based features for machine-learning (ML) studies can be engineered via the sub-packages
"featurize" and "structuregraphs". Alongside its Python interface, it also provides an easy-to-use command line
interface (CLI) that runs automatic analysis of the computations and generates a summary of results and publication-ready figures.

_LobsterPy_ has been used to produce the results in [@ngo2023dft; @morgan2023structures; @naik2023quantumchemical] and is also part of
[@atomate2] bonding analysis workflow for generating bonding analysis summaries.

# Statement of need
Although notion of "bonds" might seem unusual from a physicist's point of view, chemists have been employing it routinely to
explain various chemical phenomena and material properties.[@hoffmann1987chemistry; @burdett1995chemical; @das2023strong; @ertural2022first; @hu2023mechanism; @dronskowski2023chemical] With the recent advances in
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

# Availability
_LobsterPy_ can also be found on [PyPI](https://pypi.org/project/lobsterpy/). Detailed software documentation,
including [installation instructions](https://jageo.github.io/LobsterPy/installation/index.html) and
[implementation details](https://jageo.github.io/LobsterPy/fundamentals/index.html) are provided. The package
also includes [tutorials](https://jageo.github.io/LobsterPy/tutorial/index.html) illustrating all the basic and advanced functionalities.

# Acknowledgements
The authors would like to acknowledge the Gauss Centre for Super
computing e.V. (www.gauss-centre.eu) for funding this project by
providing generous computing time on the GCS Supercomputer
SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de)
(project pn73da) that enabled rigorous testing of this
package on a diverse set of compounds. The authors thank Jonas Grandel for reviewing the docstrings and testing package functionalities and tutorials. The authors would also like to acknowledge the maintainers of pymatgen and LOBSTER program developers.

# References
