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
Lobsterpy is a Python package developed to systematically analyze, 
describe, and visualize LOBSTER computations results. It also provides 
an easy-to-use command line interface that runs automatic analysis of 
the computations and generates a summary of results and publication-ready 
figures alongside its Python interface. Since its first release, its 
capabilities have been extended significantly. Lobsterpy can now also 
analyze crystal orbital overlap populations (COOP)[@hughbanks1983chains] and Crystal orbital 
bond index (COBI)[@müller2021crystal] files to extract summarized bonding information 
comprising coordination environments bond strengths, most relevant bonds, 
and their corresponding bonding and anti-bonding contributions. Obtaining
the most relevant orbital interaction information in each detected bond is
also possible. It also includes featurize and structure graphs modules that 
provide a pathway to get features to be used further for machine learning studies. 
The features section comprehensively overviews previously existing and newly 
added functionalities. 

Lobsterpy was used to produce the results in [@naik2023quantumchemical]

# Statement of need
Although the idea of "chemical bonding" might seem perplexing from a 
physical standpoint, it has been employed several times to explain 
various chemical phenomena and material properties. With the recent 
advances in automation frameworks for high-throughput computational 
investigations, bonding analysis for thousands of crystalline materials 
could be performed with few lines of code.[@george2022automated] This 
automation helps reduce the common mistakes inexperienced users make 
while performing bonding analysis. However, to make the inputs consistent 
and interpret the results correctly, it is necessary to systematically 
generate inputs and analyze the vast number of output files. Furthermore, 
having data from high-throughput calculations ready to utilize as inputs 
would benefit data-driven material science research.

# Features

# Availability
Lobsterpy can be found on GitHub and is also available from PyPi. 
Software documentation and installation instructions, along with 
tutorials, can be found in the documentation. The package also 
comes with several Jupyter Notebook tutorials illustrating the 
features. Several examples for demonstrating the cli features 
are also included and can be found in the software documentation.


# Acknowledgements
The authors would like to acknowledge the Gauss Centre for Super 
computing e.V. (www.gauss-centre.eu) for funding this project by 
providing generous computing time on the GCS Supercomputer 
SuperMUC-NG at Leibniz Supercomputing Centre (www.lrz.de) 
(project pn73da) that enabled rigorous testing of this 
package on a diverse set of compounds. We also acknowledge 
the maintainers of pymatgen.

# References
