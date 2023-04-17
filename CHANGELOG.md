## v0.2.9
- fix the error handling in cases ChemEnv cannot determine a coordination environment and we use coordination numbers instead for the cation-anion mode

## v0.2.8
- fix issue while saving files

## v0.2.7
- compatible with atomate2
## v0.2.5
- compatible with latest pymatgen version
## v0.2.4
- fixing linting errors due to new mypy versions
## v0.2.3
- bonding and antibonding contributions will now be integrated and a percentage of antibonding interactions below Efermi will be given.

## v0.2.2
- users can provide their own basis functions for lobsterin/INCAR generation
- documentation added
- fixes to saving files

## v0.2.1
- Fix error message when LobsterPy is used in cation-anion mode for materials that are not ionic.
- automatic plots are now saved correctly. Before only the last plot was saved.
- Additional Gaussian broadening available for COHPs
- lobsterins can be generated with the command-line interface
## v0.1.0
- First LobsterPy release
- Automatic COHP analysis (description and plots)
- Command line tool to perform automatic analysis
- Command line tool to plot COHPs, COOPs, and COBIs - also orbitalwise and summed
- Many options to refine plots (own matplotlib styles, changes of font size, sizes)
