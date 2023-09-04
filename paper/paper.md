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
Although the idea of "chemical bonding" might seem perplexing from a physical standpoint, it has been employed on several occasions to explain various chemical phenomena and material properties.[@naik2023quantumchemical]  With the recent advances in automation frameworks for high-throughput computational investigations, bonding analysis for thousands of crystalline materials could be performed with few lines of code.[@george2022automated] This automation helps to reduce the common mistakes made by inexperienced users while doing bonding analysis. However, in order to understand the results properly, it is necessary to analyze the large number of output files systematically.  Furthermore, having data from high-throughput calculations ready to utilize as inputs would benefit data-driven material science research.

`Lobsterpy` is a Python package that aims to address this need and can be used to analyze, describe, and visualize LOBSTER computations results systematically. It includes an optional featurize module that engineers feature from the output files that could be used for machine learning studies. It also has an easy-to-use command line interface that runs automatic analysis of the computations and generates a summary of results and publication-ready figures.

# Implementation and Features

# Acknowledgements

# References
