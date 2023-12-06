# Changelog

## v0.3.4
- fix changelog by @JaGeo in https://github.com/JaGeo/LobsterPy/pull/188
- Update documentation by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/180

**Full Changelog**: https://github.com/JaGeo/LobsterPy/compare/v0.3.3...v0.3.4

## v0.3.3
- fixing which_bonds  by @JonasGrandel in https://github.com/JaGeo/LobsterPy/pull/168
- fix create inputs alias not working; update test for the same by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/171
- Automatic orbital wise analysis functionality in analyze module by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/132
- Replace unittests with pytests + update CI workflow and code doc-strings by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/177
- replace get_anion_types with anion_types by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/174
- rename keys of calc quality summary and snake_case by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/178
- Remove read the docs config  by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/182
- Refactor cli.py for cleaner options on cli help by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/172
- overwrite add_cohp for interactive plotter > Now it works as expected by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/176
- Extend featurizer by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/179
- Update README.md by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/184
- Increase test coverage by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/185
- add POSCAR.lobster support in featurizer by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/186

**Full Changelog**: https://github.com/JaGeo/LobsterPy/compare/v0.3.2...v0.3.3

## v0.3.2
- cli invert axis, add get site all orbitals dos plot by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/149
- Update README.md by @JaGeo in https://github.com/JaGeo/LobsterPy/pull/150
- Bump mendeleev from 0.12.1 to 0.14.0 by @dependabot in https://github.com/JaGeo/LobsterPy/pull/151
- snakecase key names of calc quality summary dict by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/157
- Fix missing matplotlib style file in package installation  by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/159
- snakecase whichbonds arg and update tests by @naik-aakash in https://github.com/JaGeo/LobsterPy/pull/161

**Full Changelog**: https://github.com/JaGeo/LobsterPy/compare/v0.3.1...v0.3.2

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
