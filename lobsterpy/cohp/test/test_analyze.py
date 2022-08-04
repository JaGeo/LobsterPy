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

        self.analyse_NaCl_comp_range = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/NaCl_comp_range/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_nan = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
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

        self.analyse_NaSbF6 = Analysis(
            path_to_poscar=TestDir / "TestData/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/NaSbF6/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaSbF6_anbd = Analysis(
            path_to_poscar=TestDir / "TestData/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/NaSbF6/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-5.5,
        )

        self.analyse_CdF = Analysis(
            path_to_poscar=TestDir / "TestData/CdF/POSCAR",
            path_to_cohpcar=TestDir / "TestData/CdF/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/CdF/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/CdF/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
        )

        self.analyse_CdF_comp_range = Analysis(
            path_to_poscar=TestDir / "TestData/CdF_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/CdF_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/CdF_comp_range/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
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
        with self.assertRaises(ValueError) as err:
            self.analyse_C = Analysis(
                path_to_poscar=TestDir / "TestData/C/POSCAR",
                path_to_cohpcar=TestDir / "TestData/C/COHPCAR.lobster",
                path_to_icohplist=TestDir / "TestData/C/ICOHPLIST.lobster",
                path_to_charge=TestDir / "TestData/C/CHARGE.lobster",
                whichbonds="cation-anion",
                cutoff_icohp=0.1,
            )
        self.assertEqual(
            err.exception.__str__(),
            "Consider switching to an analysis of all bonds and not only cation-anion bonds. It looks like no cations are detected.",
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

    def test_all_attributes_NaSbF6(self):
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["formula"], "NaSbF6"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.98893,
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["number_of_considered_ions"],
            2,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"][
                    "F"
                ]["ICOHP_sum"]
            ),
            -3.65,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "bonding"
            ]["integral"],
            3.77,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "bonding"
            ]["perc"],
            0.95929,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "antibonding"
            ]["perc"],
            0.04071,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "antibonding"
            ]["integral"],
            0.16,
        )
        self.assertAlmostEqual(
            abs(
                float(
                    self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"][
                        "F"
                    ]["ICOHP_sum"]
                )
            ),
            (
                (
                    round(
                        self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0][
                            "bonds"
                        ]["F"]["bonding"]["integral"]
                        - self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0][
                            "bonds"
                        ]["F"]["antibonding"]["integral"],
                        2,
                    )
                )
            ),
            delta=0.10,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["ion"], "Na"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["charge"], 0.91
        )
        self.assertListEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0][
                "relevant_bonds"
            ],
            ["21", "25", "31", "34", "43", "47"],
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"][
                    "F"
                ]["ICOHP_sum"]
            ),
            -32.71,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"][
                "has_antibdg_states_below_Efermi"
            ],
            False,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"][
                "bonding"
            ]["perc"],
            1.0,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"][
                "antibonding"
            ]["perc"],
            0.0,
        )
        self.assertAlmostEqual(
            abs(
                float(
                    self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"][
                        "F"
                    ]["ICOHP_sum"]
                )
            ),
            (
                (
                    round(
                        self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1][
                            "bonds"
                        ]["F"]["bonding"]["integral"]
                        - self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1][
                            "bonds"
                        ]["F"]["antibonding"]["integral"],
                        2,
                    )
                )
            ),
            delta=0.10,
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["ion"], "Sb"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["charge"], 2.91
        )
        self.assertListEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1][
                "relevant_bonds"
            ],
            ["63", "69", "73", "80", "81", "87"],
        )
        self.assertEqual(
            self.analyse_NaSbF6.condensed_bonding_analysis["type_charges"], "Mulliken"
        )

    def test_all_attributes_NaSbF6_anbd(self):
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["formula"], "NaSbF6"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.98893,
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis[
                "number_of_considered_ions"
            ],
            2,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["env"],
            "O:6",
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0][
                    "bonds"
                ]["F"]["ICOHP_sum"]
            ),
            -3.65,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"][
                "F"
            ]["has_antibdg_states_below_Efermi"],
            True,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"][
                "F"
            ]["number_of_bonds"],
            6,
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"][
                "F"
            ]["bonding"]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"][
                "F"
            ]["antibonding"]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["ion"], "Na"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["charge"],
            0.91,
        )
        self.assertListEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0][
                "relevant_bonds"
            ],
            ["21", "25", "31", "34", "43", "47"],
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["env"],
            "O:6",
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1][
                    "bonds"
                ]["F"]["ICOHP_sum"]
            ),
            -32.71,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"][
                "F"
            ]["has_antibdg_states_below_Efermi"],
            False,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"][
                "F"
            ]["number_of_bonds"],
            6,
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"][
                "F"
            ]["bonding"]["perc"],
            1.0,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"][
                "F"
            ]["antibonding"]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["ion"], "Sb"
        )
        self.assertAlmostEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["charge"],
            2.91,
        )
        self.assertListEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1][
                "relevant_bonds"
            ],
            ["63", "69", "73", "80", "81", "87"],
        )
        self.assertEqual(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["type_charges"],
            "Mulliken",
        )

    def test_all_attributes_CdF(self):
        self.assertEqual(self.analyse_CdF.condensed_bonding_analysis["formula"], "CdF2")
        self.assertAlmostEqual(
            self.analyse_CdF.condensed_bonding_analysis["max_considered_bond_length"],
            5.99501,
        )
        self.assertAlmostEqual(
            self.analyse_CdF.condensed_bonding_analysis["number_of_considered_ions"],
            1,
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["env"], "C:8"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                    "ICOHP_sum"
                ]
            ),
            -4.96,
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "number_of_bonds"
            ],
            8,
        )
        self.assertAlmostEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "bonding"
            ]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
                "antibonding"
            ]["perc"],
            1.0,
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["ion"], "Cd"
        )
        self.assertAlmostEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["charge"], 1.57
        )
        self.assertListEqual(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["relevant_bonds"],
            ["25", "32", "35", "36", "57", "58", "61", "68"],
        )
        self.assertEqual(
            self.analyse_CdF.condensed_bonding_analysis["type_charges"], "Mulliken"
        )

    def test_all_attributes_NaCl_nan(self):
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["formula"], "NaCl"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis[
                "max_considered_bond_length"
            ],
            5.69169,
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis[
                "number_of_considered_ions"
            ],
            1,
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["env"], "O:6"
        )
        self.assertAlmostEqual(
            float(
                self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"][
                    "Cl"
                ]["ICOHP_sum"]
            ),
            -3.39,
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "has_antibdg_states_below_Efermi"
            ],
            True,
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "number_of_bonds"
            ],
            6,
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "bonding"
            ]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "antibonding"
            ]["perc"],
            0.0,
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["ion"], "Na"
        )
        self.assertAlmostEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["charge"], 0.78
        )
        self.assertListEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0][
                "relevant_bonds"
            ],
            ["21", "23", "24", "27", "28", "30"],
        )
        self.assertEqual(
            self.analyse_NaCl_nan.condensed_bonding_analysis["type_charges"], "Mulliken"
        )

    def test_ICOHP_sum_NaCl(self):
        self.assertAlmostEqual(
            abs(
                float(
                    self.analyse_NaCl_comp_range.condensed_bonding_analysis["sites"][0][
                        "bonds"
                    ]["Cl"]["ICOHP_sum"]
                )
            ),
            (
                (
                    round(
                        self.analyse_NaCl_comp_range.condensed_bonding_analysis[
                            "sites"
                        ][0]["bonds"]["Cl"]["bonding"]["integral"]
                        - self.analyse_NaCl_comp_range.condensed_bonding_analysis[
                            "sites"
                        ][0]["bonds"]["Cl"]["antibonding"]["integral"],
                        2,
                    )
                )
            ),
            delta=0.10,
        )

    def test_ICOHP_sum_CdF(self):
        self.assertAlmostEqual(
            abs(
                float(
                    self.analyse_CdF_comp_range.condensed_bonding_analysis["sites"][0][
                        "bonds"
                    ]["F"]["ICOHP_sum"]
                )
            ),
            (
                (
                    round(
                        self.analyse_CdF_comp_range.condensed_bonding_analysis["sites"][
                            0
                        ]["bonds"]["F"]["bonding"]["integral"]
                        - self.analyse_CdF_comp_range.condensed_bonding_analysis[
                            "sites"
                        ][0]["bonds"]["F"]["antibonding"]["integral"],
                        2,
                    )
                )
            ),
            delta=0.10,
        )


if __name__ == "__main__":
    unittest.main()
