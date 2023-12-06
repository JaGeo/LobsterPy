description-quality
===================

Deliver a text description of the LOBSTER calc quality analysis. Mandatory required files: POSCAR, POTCAR or POTCAR symbols, lobsterout, lobsterin. Optional files (BVA comparison): CHARGE.lobster, (DOS comparison): DOSCAR.lobster/ DOSCAR.LSO.lobster, Vasprun.xml.

.. argparse::
   :module: lobsterpy.cli
   :func: get_parser
   :prog: lobsterpy
   :path: description-quality
