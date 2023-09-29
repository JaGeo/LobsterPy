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
            bonds="all",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            n_jobs=3,
        )

        self.summary_featurize_without_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            n_jobs=3,
        )

        self.summary_featurize_with_json_overall = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="overall",
            include_cobi_data=True,
            include_coop_data=True,
            e_range=[-15, 0],
            n_jobs=3,
        )

        self.summary_featurize_with_json_bonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="bonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="mulliken",
            n_jobs=3,
        )

        self.summary_featurize_with_json_antibonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            bonds="cation_anion",
            path_to_jsons=TestDir / "TestData/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="loewdin",
            n_jobs=3,
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

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

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

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

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

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

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
            n_jobs=3,
        )

        self.fp_cohp_bonding = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=False,
            tanimoto=True,
            n_jobs=3,
        )

        self.fp_cobi = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="antibonding",
            normalize=True,
            tanimoto=True,
            fingerprint_for="cobi",
            n_jobs=3,
        )

        self.fp_coop = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "TestData/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=True,
            tanimoto=False,
            fingerprint_for="coop",
            n_jobs=3,
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

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0.000017, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0.000000, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        self.assertAlmostEqual(df.loc["mp-1000", "mp-2176"], 0.001532, places=5)

    def test_fp_cobi(self):
        fp_df = self.fp_cobi.fingerprint_df
        df = self.fp_cobi.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        self.assertAlmostEqual(df.loc["mp-1000", "mp-2176"], 0, places=5)

    def test_fp_coop(self):
        fp_df = self.fp_coop.fingerprint_df
        df = self.fp_coop.get_similarity_matrix_df()

        self.assertAlmostEqual(df.loc["mp-463", "mp-1000"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-2176"], 0, places=5)
        self.assertAlmostEqual(df.loc["mp-463", "mp-463"], 1, places=5)
        self.assertAlmostEqual(df.loc["mp-1000", "mp-2176"], 0, places=5)


class TestExceptions(unittest.TestCase):
    def test_batch_summary_featurizer_exception(self):
        with self.assertRaises(Exception) as err1:
            self.summary_featurize_with_json = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                / "TestData/Featurizer_test_data/Lobster_calcs_exceptions/1/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json.get_df()

        self.assertEqual(
            err1.exception.__str__(),
            "COBICAR.lobster or ICOBILIST.lobster file not found in mp-2176",
        )

        with self.assertRaises(Exception) as err2:
            self.summary_featurize_with_json = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                / "TestData/Featurizer_test_data/Lobster_calcs_exceptions/2/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json.get_df()

        self.assertEqual(
            err2.exception.__str__(),
            "COOPCAR.lobster or ICOOPLIST.lobster file not found in mp-1000",
        )

        # COXX exception
        with self.assertRaises(Exception) as err3:
            self.raise_coxx_exception = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "TestData/JSONS/"
            )

            _ = self.raise_coxx_exception._featurizecoxx(
                path_to_lobster_calc=self.raise_coxx_exception.path_to_lobster_calcs
            )

        self.assertEqual(
            err3.exception.__str__(),
            "COHPCAR.lobster or POSCAR or ICOHPLIST.lobster file not found in JSONS",
        )

        # Charges exception
        with self.assertRaises(Exception) as err4:
            self.raise_ch_exception = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "TestData/JSONS/"
            )

            _ = self.raise_ch_exception._featurizecharges(
                path_to_lobster_calc=self.raise_ch_exception.path_to_lobster_calcs
            )

        self.assertEqual(
            err4.exception.__str__(),
            "CHARGE.lobster or POSCAR not found in JSONS",
        )

        # Fingerprint similarity exception
        with self.assertRaises(Exception) as err8:
            fp_cohp_bonding = BatchCoxxFingerprint(
                path_to_lobster_calcs=TestDir
                / "TestData/Featurizer_test_data/Lobster_calcs",
                e_range=[-15, 0],
                feature_type="bonding",
                normalize=True,
                tanimoto=True,
                n_jobs=3,
            )

            fp_df = fp_cohp_bonding.fingerprint_df

            _ = fp_cohp_bonding._get_fp_similarity(
                fp_df.loc["mp-1000", "COXX_FP"],
                fp_df.loc["mp-2176", "COXX_FP"],
                tanimoto=True,
                normalize=True,
            )

        self.assertEqual(
            err8.exception.__str__(),
            "Cannot compute similarity index. Please set either normalize=True or "
            "tanimoto=True or both to False.",
        )


if __name__ == "__main__":
    unittest.main()
