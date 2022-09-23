Python interface
================

Import the necessary modules

.. code:: ipython3

    import os
    from lobsterpy.cohp.analyze import Analysis
    from lobsterpy.cohp.describe import Description
    import warnings
    warnings.filterwarnings('ignore')

.. code:: ipython3

    directory = "Basis_0" #Directory of your VASP and Lobster computations
    # Setup analysis dict
    analyse = Analysis(
        path_to_poscar=os.path.join(directory, "POSCAR"),
        path_to_icohplist=os.path.join(directory, "ICOHPLIST.lobster"),
        path_to_cohpcar=os.path.join(directory, "COHPCAR.lobster"),
        path_to_charge=os.path.join(directory, "CHARGE.lobster"),
        whichbonds="cation-anion",
    )

.. code:: ipython3

    # Setup Desciption dict
    describe = Description(analysis_object=analyse)
    describe.write_description()


.. parsed-literal::

    The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.
    Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOHP: -0.62 eV, antibonding interaction below EFermi) bonds.


.. code:: ipython3

    # Automatic plots
    describe.plot_cohps(ylim=[-10, 2], xlim=[-4, 4])



.. image:: Lobsterpy_tutorial_files/Lobsterpy_tutorial_38_0.png


.. code:: ipython3

    print(analyse.condensed_bonding_analysis) # dicts that summarize the results


.. parsed-literal::

    {'formula': 'CdF2', 'max_considered_bond_length': 5.98538, 'limit_icohp': (-100000, -0.1), 'number_of_considered_ions': 1, 'sites': {0: {'env': 'C:8', 'bonds': {'F': {'ICOHP_mean': '-0.62', 'ICOHP_sum': '-4.97', 'has_antibdg_states_below_Efermi': True, 'number_of_bonds': 8}}, 'ion': 'Cd', 'charge': 1.57, 'relevant_bonds': ['29', '30', '33', '40', '53', '60', '63', '64']}}, 'type_charges': 'Mulliken'}


.. code:: ipython3

    print(analyse.final_dict_bonds) # dicts that summarize the results


.. parsed-literal::

    {'Cd-F': {'ICOHP_mean': -0.62125, 'has_antbdg': True}}


.. code:: ipython3

    print(analyse.final_dict_ions) # dicts that summarize the results


.. parsed-literal::

    {'Cd': {'C:8': 1}}


