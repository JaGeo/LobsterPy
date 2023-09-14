## v0.3.1
- ICOHP vs bond length plotter by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/137
- Remove python 3.8 support by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/134
- Add units to plotters by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/136
- Bump sphinx-pdj-theme from 0.2.1 to 0.4.0 by @dependabot in https://github.com/JaGeo/LobsterPy/pull/111
- Doscar plotting by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/138
- added lobster calc quality summary method to analyze module by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/115
- Add featurizer by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/102
- COBI COOP extension by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/128
- Structure graphs by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/63
- fix docs build, remove unwanted dependencies by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/146

**Full Changelog**: https://github.com/JaGeo/LobsterPy/compare/v0.3.0...v0.3.1

## v0.3.0
- addition of an interactive plotter by @naik-aakash and @kaueltzen. Reviews by @ajjackson and @jageo
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
