Tutorial
=========

Written by Aakash Naik (aakash.naik@bam.de) based on suggestions provided by Prof. Dr. Janine George (janine.george@bam.de)


Prerequisites

* Familiar with `LOBSTER software <http://www.cohp.de>`_ and its output files
* Familiar with the :doc:`../fundamentals/index` guide. If not, please read that guide before proceeding.
* Basic Python knowledge


This tutorial will demonstrate how to use the LobsterPy package using example code snippets. By the end of this tutorial, you will be able to use the package to:

1. Automatically analyze the lobster outputs using the Python interface
2. Generate custom plots via plotting utilities 
3. Use structuregraph to generate graph objects consisting of the LOBSTER data
4. Use featurizers to extract the LOBSTER bonding analysis data as features for the ML studies
5. Get automatic analysis results and plots using command line utilities

**Note**: The following tutorial is jupyter Notebook, you can download it and execute locally. Please ignore the code blocks which have `remove-cell` tags. You can see these cell tags using `View -> Cell Toolbar -> Tags` (These code blocks are hidden in rendered docs to keep it consistent with example files paths you will use.)

.. toctree::
   :maxdepth: 2

   tutorial
   commandlineinterface
   atomateauto
   computingtimes
