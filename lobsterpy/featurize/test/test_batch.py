import os
import sys
import unittest
import pandas as pd
from pathlib import Path
from lobsterpy.featurize.batch import BatchSummaryFeaturizer, BatchCoxxFingerprint

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


class TestBatchSummaryFeaturizer(unittest.TestCase):
    def setUp(self):
        self.summary_featurize_with_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all_bonds",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
        )

        self.summary_featurize_without_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all_bonds",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
        )

        self.summary_featurize_with_json_overall = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all_bonds",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="overall",
            include_cobi_data=True,
            include_coop_data=True,
            e_range=[-15, 0],
        )

        self.summary_featurize_with_json_bonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all_bonds",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="bonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="mulliken",
        )

        self.summary_featurize_with_json_antibonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="cation_anion_bonds",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="loewdin",
        )

    def test_summary_featurize_with_json(self):
        df = self.summary_featurize_with_json.get_df()

        self.assertIsInstance(df, pd.DataFrame)

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
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        self.assertEqual(list(df.columns), expected_cols)

        expected_index = ["mp-463", "mp-1000", "mp-2176"]

        self.assertEqual(list(df.index), expected_index)

    def test_summary_featurize_without_json(self):
        df = self.summary_featurize_without_json.get_df()

        self.assertIsInstance(df, pd.DataFrame)

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
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        self.assertEqual(list(df.columns), expected_cols)

        expected_index = ["mp-463", "mp-1000", "mp-2176"]

        self.assertEqual(list(df.index), expected_index)

    def test_summary_featurize_with_json_overall(self):
        df = self.summary_featurize_with_json_overall.get_df()

        self.assertIsInstance(df, pd.DataFrame)

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
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "bnd_wICOBI",
            "antibnd_wICOBI",
            "w_ICOBI",
            "EIN_ICOBI",
            "center_COBI",
            "width_COBI",
            "skewness_COBI",
            "kurtosis_COBI",
            "bnd_wICOOP",
            "antibnd_wICOOP",
            "w_ICOOP",
            "EIN_ICOOP",
            "center_COOP",
            "width_COOP",
            "skewness_COOP",
            "kurtosis_COOP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        self.assertEqual(list(df.columns), expected_cols)

        expected_index = ["mp-463", "mp-1000", "mp-2176"]

        self.assertEqual(list(df.index), expected_index)

    def test_summary_featurize_with_json_bonding(self):
        df = self.summary_featurize_with_json_bonding.get_df()

        self.assertIsInstance(df, pd.DataFrame)

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
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "Ionicity_Mull",
        ]

        self.assertEqual(list(df.columns), expected_cols)

    def test_summary_featurize_with_json_antibonding(self):
        df = self.summary_featurize_with_json_antibonding.get_df()

        self.assertIsInstance(df, pd.DataFrame)

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
            "bnd_wICOHP",
            "antibnd_wICOHP",
            "w_ICOHP",
            "EIN_ICOHP",
            "center_COHP",
            "width_COHP",
            "skewness_COHP",
            "kurtosis_COHP",
            "Ionicity_Loew",
        ]

        self.assertEqual(list(df.columns), expected_cols)


class TestBatchCoxxFingerprint(unittest.TestCase):
    def setUp(self):
        self.fp_cohp_overall = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="overall",
            normalize=True,
            tanimoto=True,
        )

        self.fp_cohp_bonding = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=False,
            tanimoto=True,
        )

        self.fp_cobi = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="antibonding",
            normalize=True,
            tanimoto=True,
            are_cobis=True,
        )

        self.fp_coop = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=True,
            tanimoto=False,
            are_coops=True,
        )

    def test_fp_cohp_overall(self):
        df = self.fp_cohp_overall.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], -0.033251, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], -0.013751, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        self.assertAlmostEqual(df.loc["mp-1000", "mp-2176"], 0.046889, places=5)

    def test_fp_cohp_bonding(self):
        fp_df = self.fp_cohp_bonding.fingerprint_df
        df = self.fp_cohp_bonding.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0.0000171, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0.000000, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        for val in fp_df.loc["mp-1000", "COXX_FP"].coxx:
            self.assertLessEqual(val, 0)

    def test_fp_cobi(self):
        fp_df = self.fp_cobi.fingerprint_df
        df = self.fp_cobi.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        for val in fp_df.loc["mp-463", "COXX_FP"].coxx:
            self.assertLessEqual(val, 0)

    def test_fp_coop(self):
        fp_df = self.fp_coop.fingerprint_df
        df = self.fp_coop.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        for val in fp_df.loc["mp-2176", "COXX_FP"].coxx:
            self.assertGreaterEqual(val, 0)


class TestExceptions(unittest.TestCase):
    def test_batch_summary_featurizer_exception(self):
        with self.assertRaises(Exception) as err:
            self.summary_featurize_with_json = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                                      / "TestData/Featurizer_test_data/Lobster_calcs_exceptions/1/",
                bonds="all_bonds",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            df = self.summary_featurize_with_json.get_df()
            print(df.to_string())
        self.assertEqual(
            err.exception.__str__(),
            "COBICAR.lobster.gz or ICOBILIST.lobster.gz file not found in mp-2176",
        )

        with self.assertRaises(Exception) as err:
            self.summary_featurize_with_json = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                                      / "TestData/Featurizer_test_data/Lobster_calcs_exceptions/2/",
                bonds="all_bonds",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            df = self.summary_featurize_with_json.get_df()
            print(df.to_string())
        self.assertEqual(
            err.exception.__str__(),
            "COOPCAR.lobster.gz or ICOOPLIST.lobster.gz file not found in mp-1000",
        )


if __name__ == "__main__":
    unittest.main()
