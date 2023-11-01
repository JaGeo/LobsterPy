import unittest
import pandas as pd
import numpy as np
from pathlib import Path
from lobsterpy.featurize.core import FeaturizeLobsterpy, FeaturizeCharges, FeaturizeCOXX

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


class TestFeaturizeLobsterpy(unittest.TestCase):
    def setUp(self):
        self.featurize_mp1249_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "TestData/JSONS/mp-1249.json.gz", bonds="all"
        )

        self.featurize_mp1249_json_ca = FeaturizeLobsterpy(
            path_to_json=TestDir / "TestData/JSONS/mp-1249.json.gz",
            bonds="cation-anion",
        )
        self.featurize_mp1958_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "TestData/JSONS/mp-1958.json.gz", bonds="all"
        )
        self.featurize_mp14652_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "TestData/JSONS/mp-14652.json.gz", bonds="all"
        )

        self.featurize_CsH_madelung = FeaturizeLobsterpy(
            path_to_lobster_calc=TestDir / "TestData/CsH/", bonds="all"
        )

    def test_featurize_mp1249_json(self):
        df = self.featurize_mp1249_json.get_df(ids="mp-1249")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Icohp_mean_avg",
            "Icohp_mean_max",
            "Icohp_mean_min",
            "Icohp_mean_std",
            "Icohp_sum_avg",
            "Icohp_sum_max",
            "Icohp_sum_min",
            "Icohp_sum_std",
            "bonding_perc_avg",
            "bonding_perc_max",
            "bonding_perc_min",
            "bonding_perc_std",
            "antibonding_perc_avg",
            "antibonding_perc_min",
            "antibonding_perc_max",
            "antibonding_perc_std",
            "Madelung_Mull",
            "Madelung_Loew",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "mp-1249")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["mp-1249", "Icohp_mean_avg"], -1.020000, places=5)
        self.assertAlmostEqual(df.loc["mp-1249", "Icohp_mean_max"], -1.020000, places=5)
        self.assertAlmostEqual(df.loc["mp-1249", "Icohp_mean_min"], -1.020000, places=5)
        self.assertAlmostEqual(df.loc["mp-1249", "Icohp_mean_std"], 0.000000, places=5)
        self.assertAlmostEqual(df.loc["mp-1249", "Madelung_Mull"], -52.000000, places=5)
        self.assertAlmostEqual(
            df.loc["mp-1249", "bonding_perc_avg"], 0.978985, places=5
        )

    def test_featurize_mp1249_json_ca(self):
        df = self.featurize_mp1249_json_ca.get_df(ids="mp-1249")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Icohp_mean_avg",
            "Icohp_mean_max",
            "Icohp_mean_min",
            "Icohp_mean_std",
            "Icohp_sum_avg",
            "Icohp_sum_max",
            "Icohp_sum_min",
            "Icohp_sum_std",
            "bonding_perc_avg",
            "bonding_perc_max",
            "bonding_perc_min",
            "bonding_perc_std",
            "antibonding_perc_avg",
            "antibonding_perc_min",
            "antibonding_perc_max",
            "antibonding_perc_std",
            "Madelung_Mull",
            "Madelung_Loew",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "mp-1249")

    def test_featurize_mp1958_json(self):
        df = self.featurize_mp1958_json.get_df()

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Icohp_mean_avg",
            "Icohp_mean_max",
            "Icohp_mean_min",
            "Icohp_mean_std",
            "Icohp_sum_avg",
            "Icohp_sum_max",
            "Icohp_sum_min",
            "Icohp_sum_std",
            "bonding_perc_avg",
            "bonding_perc_max",
            "bonding_perc_min",
            "bonding_perc_std",
            "antibonding_perc_avg",
            "antibonding_perc_min",
            "antibonding_perc_max",
            "antibonding_perc_std",
            "Madelung_Mull",
            "Madelung_Loew",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "mp-1958")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["mp-1958", "Icohp_sum_avg"], -2.96000, places=5)
        self.assertAlmostEqual(df.loc["mp-1958", "Icohp_sum_max"], -2.96000, places=5)
        self.assertAlmostEqual(df.loc["mp-1958", "Icohp_sum_min"], -2.96000, places=5)
        self.assertAlmostEqual(df.loc["mp-1958", "Icohp_sum_std"], 0.000000, places=5)
        self.assertAlmostEqual(df.loc["mp-1958", "Madelung_Loew"], -16.68000, places=5)
        self.assertAlmostEqual(
            df.loc["mp-1958", "antibonding_perc_avg"], 0.14528, places=5
        )

    def test_featurize_mp14652_json(self):
        df = self.featurize_mp14652_json.get_df(ids="mp-14652")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Icohp_mean_avg",
            "Icohp_mean_max",
            "Icohp_mean_min",
            "Icohp_mean_std",
            "Icohp_sum_avg",
            "Icohp_sum_max",
            "Icohp_sum_min",
            "Icohp_sum_std",
            "bonding_perc_avg",
            "bonding_perc_max",
            "bonding_perc_min",
            "bonding_perc_std",
            "antibonding_perc_avg",
            "antibonding_perc_min",
            "antibonding_perc_max",
            "antibonding_perc_std",
            "Madelung_Mull",
            "Madelung_Loew",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "mp-14652")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["mp-14652", "Icohp_mean_std"], 2.335070, places=5)
        self.assertAlmostEqual(
            df.loc["mp-14652", "bonding_perc_max"], 0.889620, places=5
        )
        self.assertAlmostEqual(
            df.loc["mp-14652", "bonding_perc_min"], 0.873420, places=5
        )
        self.assertAlmostEqual(
            df.loc["mp-14652", "bonding_perc_std"], 0.006339, places=5
        )
        self.assertAlmostEqual(
            df.loc["mp-14652", "antibonding_perc_min"], 0.110380, places=5
        )
        self.assertAlmostEqual(
            df.loc["mp-14652", "antibonding_perc_max"], 0.126580, places=5
        )
        self.assertAlmostEqual(
            df.loc["mp-14652", "antibonding_perc_std"], 0.006339, places=5
        )

    def test_featurize_CsH_madelung(self):
        df = self.featurize_CsH_madelung.get_df()

        self.assertTrue(np.isnan(df.loc["CsH", "Madelung_Mull"]))
        self.assertTrue(np.isnan(df.loc["CsH", "Madelung_Loew"]))


class TestFeaturizeCOXX(unittest.TestCase):
    def setUp(self):
        self.featurize_NaCl_COXX = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icoxxlist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_structure=TestDir / "TestData/NaCl/POSCAR",
            feature_type="overall",
            e_range=[-5, 0],
        )
        self.featurize_CdF_COXX = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "TestData/CdF/COHPCAR.lobster",
            path_to_icoxxlist=TestDir / "TestData/CdF/ICOHPLIST.lobster",
            path_to_structure=TestDir / "TestData/CdF/POSCAR",
            feature_type="bonding",
            e_range=[-5, 0],
        )
        self.featurize_K3Sb_COXX = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "TestData/K3Sb/COHPCAR.lobster.gz",
            path_to_icoxxlist=TestDir / "TestData/K3Sb/ICOHPLIST.lobster.gz",
            path_to_structure=TestDir / "TestData/K3Sb/POSCAR.gz",
            feature_type="antibonding",
            e_range=[-5, 0],
        )

    def test_featurize_NaCl_COXX(self):
        df = self.featurize_NaCl_COXX.get_summarized_coxx_df(ids="NaCl")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "edge_COHP",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "NaCl")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["NaCl", "bnd_wICOHP"], 97.233878, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "antibnd_wICOHP"], 2.766122, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "w_ICOHP"], -0.150558, places=5)

        self.assertAlmostEqual(df.loc["NaCl", "EIN_ICOHP"], 27.843536, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "center_COHP"], -4.96241, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "width_COHP"], 8.881784e-16, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "skewness_COHP"], 1, places=5)
        self.assertAlmostEqual(df.loc["NaCl", "kurtosis_COHP"], 1, places=5)

        # test summary features using label list
        df1 = self.featurize_NaCl_COXX.get_summarized_coxx_df(label_list=["2", "3"])
        self.assertNotEqual(
            df.loc["NaCl", "center_COHP"],
            df1.loc["NaCl", "center_COHP"],
        )

    def test_featurize_NaCl_COXX_fingerprint(self):
        df = self.featurize_NaCl_COXX.get_coxx_fingerprint_df(n_bins=20000)

        fingerprint = df.loc["NaCl", "COXX_FP"]

        self.assertNotEqual(fingerprint.n_bins, 20000)

        df1 = self.featurize_NaCl_COXX.get_coxx_fingerprint_df(binning=False)

        fingerprint = df1.loc["NaCl", "COXX_FP"]

        self.assertEqual(fingerprint.n_bins, 401)

        df2 = self.featurize_NaCl_COXX.get_coxx_fingerprint_df(label_list=["3", "5"])

        fingerprint_label = df2.loc["NaCl", "COXX_FP"]

        self.assertNotEqual(fingerprint.__str__(), fingerprint_label.__str__())

    def test_featurize_CdF_COXX_fingerprint(self):
        df = self.featurize_CdF_COXX.get_coxx_fingerprint_df(n_bins=20000)

        fingerprint = df.loc["CdF", "COXX_FP"]

        self.assertNotEqual(fingerprint.n_bins, 20000)

        df1 = self.featurize_CdF_COXX.get_coxx_fingerprint_df(binning=False)

        fingerprint = df1.loc["CdF", "COXX_FP"]

        self.assertEqual(fingerprint.n_bins, 401)

        df2 = self.featurize_CdF_COXX.get_coxx_fingerprint_df(label_list=["3", "5"])

        fingerprint_label = df2.loc["CdF", "COXX_FP"]

        self.assertNotEqual(fingerprint.__str__(), fingerprint_label.__str__())

    def test_featurize_CdF_COXX(self):
        df = self.featurize_CdF_COXX.get_summarized_coxx_df()

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "edge_COHP",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "CdF")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["CdF", "bnd_wICOHP"], 81.112732, places=5)
        self.assertAlmostEqual(df.loc["CdF", "antibnd_wICOHP"], 18.887268, places=5)
        self.assertAlmostEqual(df.loc["CdF", "w_ICOHP"], -0.198235, places=5)

        self.assertAlmostEqual(df.loc["CdF", "EIN_ICOHP"], 18.76634, places=5)
        self.assertAlmostEqual(df.loc["CdF", "center_COHP"], -4.748383, places=5)
        self.assertAlmostEqual(df.loc["CdF", "width_COHP"], 0.157761, places=5)
        self.assertAlmostEqual(df.loc["CdF", "skewness_COHP"], 0.910094, places=5)
        self.assertAlmostEqual(df.loc["CdF", "kurtosis_COHP"], 2.866611, places=5)

        # test using label list
        df1 = self.featurize_CdF_COXX.get_summarized_coxx_df(
            label_list=["2", "3", "30"]
        )
        self.assertNotEqual(
            df.loc["CdF", "center_COHP"],
            df1.loc["CdF", "center_COHP"],
        )

    def test_featurize_K3Sb_COXX(self):
        df = self.featurize_K3Sb_COXX.get_summarized_coxx_df(ids="K3Sb")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "edge_COHP",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "K3Sb")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["K3Sb", "bnd_wICOHP"], 97.019044, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "antibnd_wICOHP"], 2.980956, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "w_ICOHP"], -0.318218, places=5)

        self.assertAlmostEqual(df.loc["K3Sb", "EIN_ICOHP"], 11.597595, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "center_COHP"], -0.198211, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "width_COHP"], 0.233826, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "skewness_COHP"], -1.626643, places=5)
        self.assertAlmostEqual(df.loc["K3Sb", "kurtosis_COHP"], 3.771873, places=5)


class TestFeaturizeCharges(unittest.TestCase):
    def setUp(self):
        self.featurize_C_Charge = FeaturizeCharges(
            path_to_structure=TestDir / "TestData/C/POSCAR",
            path_to_charge=TestDir / "TestData/C/CHARGE.lobster",
            charge_type="mulliken",
        )
        self.featurize_CdF_Charge = FeaturizeCharges(
            path_to_structure=TestDir / "TestData/CdF/POSCAR",
            path_to_charge=TestDir / "TestData/CdF/CHARGE.lobster",
            charge_type="mulliken",
        )
        self.featurize_K3Sb_Charge = FeaturizeCharges(
            path_to_structure=TestDir / "TestData/K3Sb/POSCAR.gz",
            path_to_charge=TestDir / "TestData/K3Sb/CHARGE.lobster.gz",
            charge_type="loewdin",
        )

    def test_featurize_C_Charge(self):
        df = self.featurize_C_Charge.get_df(ids="C")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Mull",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "C")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["C", "Ionicity_Mull"], 0.0, places=5)

    def test_featurize_CdF_Charge(self):
        df = self.featurize_CdF_Charge.get_df(ids="CdF")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Mull",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "CdF")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["CdF", "Ionicity_Mull"], 0.788333, places=5)

    def test_featurize_K3Sb_Charge(self):
        df = self.featurize_K3Sb_Charge.get_df(ids="K3Sb")

        # Test that the function returns a pandas DataFrame
        self.assertIsInstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Loew",
        ]
        self.assertCountEqual(list(df.columns), expected_cols)

        # Test that the DataFrame has the expected index
        self.assertEqual(df.index[0], "K3Sb")

        # Test that all the values in the DataFrame
        self.assertAlmostEqual(df.loc["K3Sb", "Ionicity_Loew"], 0.563333, places=5)


class TestExceptions(unittest.TestCase):
    def test_lobsterpy_featurize_exception(self):
        with self.assertRaises(Exception) as err:
            self.featurize_mp1249_json = FeaturizeLobsterpy(
                path_to_json=None, path_to_lobster_calc=None, bonds="all"
            )

            _ = self.featurize_mp1249_json.get_df()

        self.assertEqual(
            err.exception.__str__(),
            "Please provide either path to lightweight lobster jsons or path to lobster calc",
        )

        with self.assertRaises(Exception) as err:
            self.featurize_mp1249_json = FeaturizeLobsterpy(
                path_to_json=None, path_to_lobster_calc=TestDir, bonds="all"
            )

            _ = self.featurize_mp1249_json.get_df()

        self.assertEqual(
            err.exception.__str__(),
            "Path provided for Lobster calc directory seems incorrect."
            "It does not contain COHPCAR.lobster, ICOHPLIST.lobster, POSCAR and "
            "CHARGE.lobster files needed for automatic analysis using LobsterPy",
        )

        with self.assertRaises(Exception) as err:
            self.featurize_CsH_cation_anion = FeaturizeLobsterpy(
                path_to_lobster_calc=TestDir / "TestData/CsH/", bonds="cation-anion"
            )

            _ = self.featurize_CsH_cation_anion.get_df()

        self.assertEqual(
            err.exception.__str__(),
            "No cation-anion bonds detected for CsH structure. "
            "Please switch to ´all´ bonds mode",
        )

        with self.assertRaises(Exception) as err:
            self.featurize_C_cation_anion = FeaturizeLobsterpy(
                path_to_lobster_calc=TestDir / "TestData/C/", bonds="cation-anion"
            )

            _ = self.featurize_C_cation_anion.get_df()

        self.assertEqual(
            err.exception.__str__(),
            "No cation-anion bonds detected for C structure. "
            "Please switch to ´all´ bonds mode",
        )

    def test_featurize_charges(self):
        with self.assertRaises(Exception) as err:
            self.featurize_CdF_Charge = FeaturizeCharges(
                path_to_structure=TestDir / "TestData/CdF/POSCAR",
                path_to_charge=TestDir / "TestData/CdF/CHARGE.lobster",
                charge_type="Mull",
            )

            _ = self.featurize_CdF_Charge.get_df()

        self.assertEqual(
            err.exception.__str__(),
            "Please check the requested charge_type. "
            "Possible options are `Mulliken` or `Loewdin`",
        )

    def test_featurize_coxx(self):
        with self.assertRaises(Exception) as err:
            self.featurize_COXX = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
                path_to_icoxxlist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
                path_to_structure=TestDir / "TestData/NaCl/POSCAR",
                feature_type="summed",
                e_range=[-5, 0],
            )

            _ = self.featurize_COXX.get_summarized_coxx_df()

        self.assertEqual(
            err.exception.__str__(),
            "Please recheck fp_type requested argument.Possible options are bonding/antibonding/overall",
        )

        with self.assertRaises(Exception) as err2:
            self.featurize_COXX = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
                path_to_icoxxlist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
                path_to_structure=TestDir / "TestData/NaCl/POSCAR",
                feature_type="bonding",
                e_range=[-5, 0],
            )

            _ = self.featurize_COXX.get_coxx_fingerprint_df(spin_type="-1")

        self.assertEqual(
            err2.exception.__str__(),
            "Check the spin_type argument." "Possible options are summed/up/down",
        )

        with self.assertRaises(Exception) as err3:
            self.featurize_COXX = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "TestData/NaSi/COHPCAR.lobster",
                path_to_icoxxlist=TestDir / "TestData/NaSi/ICOHPLIST.lobster",
                path_to_structure=TestDir / "TestData/NaSi/POSCAR",
                feature_type="bonding",
                e_range=[-5, 0],
                are_cobis=True,
                are_coops=True,
            )

            _ = self.featurize_COXX.get_coxx_fingerprint_df()

        self.assertEqual(
            err3.exception.__str__(),
            "You cannot have info about COOPs and COBIs in the same file.",
        )

        with self.assertRaises(Exception) as err:
            self.featurize_NaCl_COXX = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
                path_to_icoxxlist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
                path_to_structure=TestDir / "TestData/NaCl/POSCAR",
                feature_type="antibond",
                e_range=[-5, 0],
            )

            _ = self.featurize_NaCl_COXX.get_summarized_coxx_df()

        self.assertEqual(
            err.exception.__str__(),
            "Please recheck fp_type requested argument."
            "Possible options are bonding/antibonding/overall",
        )

        with self.assertRaises(Exception) as err:
            self.featurize_NaCl_COXX = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
                path_to_icoxxlist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
                path_to_structure=TestDir / "TestData/NaCl/POSCAR",
                feature_type="antibonding",
                e_range=[-5, 0],
            )

            _ = self.featurize_NaCl_COXX.get_coxx_fingerprint_df(spin_type="down")

        self.assertEqual(
            err.exception.__str__(),
            "LOBSTER calculation is non-spin polarized. "
            "Please switch spin_type to `up`",
        )


if __name__ == "__main__":
    unittest.main()
