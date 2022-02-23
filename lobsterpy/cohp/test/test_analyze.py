import unittest
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


class TestAnalyse(unittest.TestCase):
    def setUp(self):
        self.analyse_NaCl = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.analyse_NaCl_valences = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=None,
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_madelung = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            path_to_madelung=TestDir / "TestData/NaCl/MadelungEnergies.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTiO3 = Analysis(
            path_to_poscar=TestDir / "TestData/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "TestData/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/BaTiO3/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTaO2N1 = Analysis(
            path_to_poscar=TestDir / "TestData/BaTaO2N1/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/BaTaO2N1/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/BaTaO2N1/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/BaTaO2N1/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTiO3_differentcutoff = Analysis(
            path_to_poscar=TestDir / "TestData/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "TestData/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/BaTiO3/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.001,
        )

        self.analyse_NaCl_distorted = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl_distorted/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl_distorted/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl_distorted/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl_distorted/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_spin = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl_spin/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl_spin/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl_spin/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl_spin/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_madelung_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            path_to_madelung=TestDir / "TestData/NaCl/MadelungEnergies.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )
        self.analyse_NaSi_madelung_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaSi/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaSi/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaSi/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaSi/CHARGE.lobster",
            path_to_madelung=TestDir / "TestData/NaSi/MadelungEnergies.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTaO2N1_cutoff = Analysis(
            path_to_poscar=TestDir / "TestData/BaTaO2N1/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/BaTaO2N1/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/BaTaO2N1/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/BaTaO2N1/CHARGE.lobster.gz",
            whichbonds="all",
            cutoff_icohp=0.001,
        )

        # different environment than O:6

    def test_exception(self):
        with self.assertRaises(ValueError):
            self.analyse_BaTaO2N1 = Analysis(
                path_to_poscar=TestDir / "TestData/BaTaO2N1/POSCAR.gz",
                path_to_cohpcar=TestDir / "TestData/BaTaO2N1/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "TestData/BaTaO2N1/ICOHPLIST.lobster.gz",
                path_to_charge=TestDir / "TestData/BaTaO2N1/CHARGE.lobster.gz",
                whichbonds="cation-cation",
                cutoff_icohp=0.1,
            )

    def test_all_attributes_NaCl_Mulliken(self):
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["formula"], "NaCl"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl.condensed_bonding_analysis["max_considered_bond_length"],
            5.69169,
        )
        self.assertAlmostEqual(
            self.analyse_NaCl.condensed_bonding_analysis["number_of_considered_ions"], 1
        )
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                    "ICOHP_mean"
                ]
            ),
            -0.57,
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                    "ICOHP_sum"
                ]
            ),
            -3.39,
        )
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["ion"], "Na"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["charge"], 0.78
        )
        self.assertListEqual(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["relevant_bonds"],
            ["21", "23", "24", "27", "28", "30"],
        )
        self.assertEqual(
            self.analyse_NaCl.condensed_bonding_analysis["type_charges"], "Mulliken"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_madelung.condensed_bonding_analysis["madelung_energy"],
            -5.40,
        )

    def test_all_attributes_NaCl_valences(self):
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["formula"], "NaCl"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.69169,
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis[
                "number_of_considered_ions"
            ],
            1,
        )
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["env"],
            "O:6",
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0][
                    "bonds"
                ]["Cl"]["ICOHP_mean"]
            ),
            -0.57,
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0][
                    "bonds"
                ]["Cl"]["ICOHP_sum"]
            ),
            -3.39,
        )
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"][
                "Cl"
            ]["has_antibdg_states_below_Efermi"],
            True,
        )
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"][
                "Cl"
            ]["number_of_bonds"],
            6,
        )
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["ion"],
            "Na",
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["charge"],
            1,
        )
        self.assertListEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0][
                "relevant_bonds"
            ],
            ["21", "23", "24", "27", "28", "30"],
        )
        self.assertEqual(
            self.analyse_NaCl_valences.condensed_bonding_analysis["type_charges"],
            "Valences",
        )

    def test_all_attributes_BaTiO3(self):
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["formula"], "BaTiO3"
        )
        self.assertAlmostEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.94882,
        )
        self.assertAlmostEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["number_of_considered_ions"],
            1,
        )
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"][
                    "O"
                ]["ICOHP_sum"]
            ),
            -21.24,
        )
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"]["O"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"]["O"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["ion"], "Ti"
        )
        self.assertAlmostEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["charge"], 0.96
        )
        self.assertListEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1][
                "relevant_bonds"
            ],
            ["87", "88", "98", "101", "109", "114"],
        )
        self.assertEqual(
            self.analyse_BaTiO3.condensed_bonding_analysis["type_charges"], "Mulliken"
        )

    def test_all_attributes_NaCl_all(self):
        print(self.analyse_NaCl_all.condensed_bonding_analysis)
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["formula"], "NaCl"
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["formula"], "NaCl"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.69169,
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis[
                "number_of_considered_ions"
            ],
            2,
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["env"], "O:6"
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][1]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"][
                    "Cl"
                ]["ICOHP_mean"]
            ),
            -0.57,
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"][
                    "Cl"
                ]["ICOHP_sum"]
            ),
            -3.39,
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["ion"], "Na"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["charge"], 0.78
        )
        self.assertListEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0][
                "relevant_bonds"
            ],
            ["21", "23", "24", "27", "28", "30"],
        )
        self.assertEqual(
            self.analyse_NaCl_all.condensed_bonding_analysis["type_charges"], "Mulliken"
        )

    def test_final_dicts(self):
        self.assertAlmostEqual(
            self.analyse_NaSi_madelung_all.final_dict_bonds["Na-Si"]["ICOHP_mean"],
            -0.44833333333,
        )
        self.assertDictEqual(self.analyse_NaCl.final_dict_ions["Na"], {"O:6": 1})
        self.assertEqual(
            self.analyse_NaCl.final_dict_bonds["Cl-Na"]["has_antbdg"], True
        )
        self.assertAlmostEqual(
            self.analyse_NaCl.final_dict_bonds["Cl-Na"]["ICOHP_mean"],
            -0.565000000000000,
        )

        self.assertDictEqual(self.analyse_BaTiO3.final_dict_ions["Ti"], {"O:6": 1})
        self.assertEqual(
            self.analyse_BaTiO3.final_dict_bonds["O-Ti"]["has_antbdg"], True
        )
        self.assertAlmostEqual(
            self.analyse_BaTiO3.final_dict_bonds["O-Ti"]["ICOHP_mean"],
            -3.5399999999999996,
        )

        self.assertDictEqual(
            self.analyse_NaCl_distorted.final_dict_ions, {"Na": {"O:6": 4}}
        )
        self.assertEqual(
            self.analyse_NaCl_distorted.final_dict_bonds["Cl-Na"]["has_antbdg"], True
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_distorted.final_dict_bonds["Cl-Na"]["ICOHP_mean"],
            -0.5658333333333334,
        )
        self.assertEqual(
            self.analyse_NaCl_spin.final_dict_bonds["Cl-Na"]["has_antbdg"], True
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_spin.final_dict_bonds["Cl-Na"]["ICOHP_mean"],
            -0.5650000000000001,
        )
        self.assertDictEqual(self.analyse_NaCl_spin.final_dict_ions, {"Na": {"O:6": 1}})
        self.assertAlmostEqual(
            self.analyse_NaCl_all.final_dict_bonds["Cl-Na"]["ICOHP_mean"],
            -0.5650000000000001,
        )
        self.assertDictEqual(
            self.analyse_NaCl_all.final_dict_ions, {"Na": {"O:6": 1}, "Cl": {"O:6": 1}}
        )


if __name__ == "__main__":
    unittest.main()
