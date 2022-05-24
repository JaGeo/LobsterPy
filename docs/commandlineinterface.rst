Command line interface
======================

Creating input files
--------------------

.. note::
   
   Important tags in INCAR of VASP to be paid attention before
   performing lobster runs are NBANDS, NSW and ISYM. It is absolutely
   necessary that VASP static run is performed (no movements of atoms,
   NSW = 0) before running lobster program. LOBSTER can only deal with
   VASP WAVECAR that contain results for the entire mesh or only half of
   it. To do this, in the INCAR set ISYM = -1 (entire mesh / symmetry
   switched off) or ISYM = 0 (half mesh/time-reversal). And to make sure
   WAVECAR are written set LWAVE = .TRUE. For pCOHP analyses one needs
   to have as many bands as there are orbitals in local basis. For pCOHP
   analyses using LOBSTER, however, you need to manually set NBANDS in
   the INCAR file.

With lobsterpy these intricate details are taken care of with single
command. We need the standard VASP input files i.e
``INCAR, KPOINTS, POTCAR and POSCAR`` in the calculation directory. Once
you have these files, one simply needs to run the following command :

``lobsterpy create-inputs``

The above command will create set of input files (INCAR and lobsterin)
depending on the basis sets that are available in Lobster.

The NBANDS, NSW, ISYM tag will be changed in existing INCAR file and new
INCAR files will be written in the existing directory. The newly created
INCAR file will be named ``INCAR.lobsterpy``\ by default. Simultaneously
``lobsterin.lobsterpy`` files are created that is necessary for lobster
run (this is the file that instructs lobster program what computations
needed to be performed).

You can also change the names of output files and path where they are
saved using following optional tags:

``lobsterpy create-inputs --incar-out <path/to/incar>/INCAR --lobsterin-out <path/to/lobsterin>/lobsterin``

You can also use help to know addtional options using
``lobsterpy create-inputs -h``

In our example ``Cd`` element has two basis sets ``4d 5s`` ``4d 5s 5p``,
thus following files are created:

::

   INCAR.lobsterpy-0
   INCAR.lobsterpy-1
   lobsterin.lobsterpy-0
   lobsterin.lobsterpy-1

The suffix “-0” & “-1” indicate input files corresponding to smaller and
larger basis of ``Cd`` respectively.

.. warning::
     
     The ‘KPOINTS’ file is not adapted, it is important for user
     to select appropriate grid density before starting VASP
     computations. Usually a factor of 50 x reciprocal lattice vectors
     is sufficient to get reliable bonding analysis results.
