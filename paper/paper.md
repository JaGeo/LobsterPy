---
title: 'LobsterPy: Package to automatically analyze Lobster runs'
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
The LOBSTER software helps extract chemical bonding information from materials by post-processing the density functional theory
computational data by projecting the wave functions onto an atomic orbital basis.  _Lobsterpy_ is a Python package that provides
convenient tools to systematically analyze, describe, and visualize these LOBSTER computations results. Since its first release,
its capabilities have been extended significantly. Unlike earlier versions which could only automatically analyze Crystal Orbital
Hamilton Populations (COHPs)[@dronskowski1993crystal], LobsterPy can now also analyze crystal orbital overlap populations
(COOP)[@hughbanks1983chains] and Crystal orbital bond index (COBI)[@müller2021crystal] to extract summarized bonding
information (includes information on coordination environments, bond strengths, most relevant bonds, bonding and anti-bonding
contributions). Optionally, users can further extract the most important orbitals contributing to the relevant bonds. Additionally,
featurize and structure graphs utility sub-packages provide a pathway to engineer features to be used further for machine learning
studies. Alongside its Python interface, it also provides an easy-to-use command line interface (CLI) that runs automatic analysis
of the computations and generates a summary of results and publication-ready figures.

_Lobsterpy_ was used to produce the results in [@ngo2023dft, @morgan2023structures, @naik2023quantumchemical]

# Statement of need
Although the idea of "chemical bonding" might seem perplexing from a
physical standpoint, it has been employed several times to explain
various chemical phenomena and material properties.[@das2023strong, @ertural2022first,
@hu2023mechanism] With the recent
advances in automation frameworks for high-throughput computational
investigations, bonding analysis for thousands of crystalline materials
could be performed with few lines of code.[@george2022automated] This
automation helps reduce the common mistakes inexperienced users make
while performing bonding analysis. However, it is essential to systematically
generate inputs and post-process the output files consistently to have
reliable and reproducible results. Furthermore,
having data from high-throughput calculations ready to utilize as inputs
would benefit data-driven material science research. _Lobsterpy_ aims to fulfill
this missing link.

# Features
- Automatic summarized bonding analysis JSONs and text descriptions based on COHPs (ICOHPs), COBIs (ICOBIs) and COOPs (ICOOPs)
- JSONs and textual description of LOBSTER calculation quality
- Static and interactive plots of most relevant COHPs, COBIs and COOPs
- Generate inputs for bonding analysis calculations
- Generate features to be used for ML studies
- Command line interface for automatic bonding analysis and plotting.


# Availability
Lobsterpy can be found on GitHub and is also available from PyPI.
Detailed software documentation including [implementation details](https://jageo.github.io/LobsterPy/fundamentals/index.html) and [installation instructions](https://jageo.github.io/LobsterPy/installation/index.html) are provided.
The package also comes with [tutorials](https://jageo.github.io/LobsterPy/tutorial/index.html) illustrating the usage and features.

# Acknowledgements
The authors would like to acknowledge the Gauss Centre for Super
computing e.V. (www.gauss-centre.eu) for funding this project by
providing generous computing time on the GCS Supercomputer
SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de)
(project pn73da) that enabled rigorous testing of this
package on a diverse set of compounds. We also acknowledge
the maintainers of pymatgen.

# References
