from pathlib import Path

import pandas as pd
import pytest

from lobsterpy.featurize.batch import BatchCoxxFingerprint, BatchSummaryFeaturizer

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestBatchSummaryFeaturizer:
    def test_summary_featurize_with_json(self):
        summary_featurize_with_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            path_to_jsons=TestDir / "test_data/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            n_jobs=3,
        )
        df = summary_featurize_with_json.get_df()

        assert isinstance(df, pd.DataFrame)

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
            "edge_COHP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        assert list(df.columns) == expected_cols

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert list(df.index) == expected_index

    def test_summary_featurize_without_json(self):
        summary_featurize_without_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            n_jobs=3,
        )

        df = summary_featurize_without_json.get_df()

        assert isinstance(df, pd.DataFrame)

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
            "edge_COHP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        assert list(df.columns) == expected_cols

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert list(df.index) == expected_index

    def test_summary_featurize_with_json_overall(self):
        summary_featurize_with_json_overall = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            path_to_jsons=TestDir / "test_data/Featurizer_test_data/JSONS",
            feature_type="overall",
            include_cobi_data=True,
            include_coop_data=True,
            e_range=[-15, 0],
            n_jobs=3,
        )

        df = summary_featurize_with_json_overall.get_df()

        assert isinstance(df, pd.DataFrame)

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
            "edge_COHP",
            "bnd_wICOBI",
            "antibnd_wICOBI",
            "w_ICOBI",
            "EIN_ICOBI",
            "center_COBI",
            "width_COBI",
            "skewness_COBI",
            "kurtosis_COBI",
            "edge_COBI",
            "bnd_wICOOP",
            "antibnd_wICOOP",
            "w_ICOOP",
            "EIN_ICOOP",
            "center_COOP",
            "width_COOP",
            "skewness_COOP",
            "kurtosis_COOP",
            "edge_COOP",
            "Ionicity_Mull",
            "Ionicity_Loew",
        ]

        assert list(df.columns) == expected_cols

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert list(df.index) == expected_index

    def test_summary_featurize_with_json_bonding(self):
        summary_featurize_with_json_bonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            path_to_jsons=TestDir / "test_data/Featurizer_test_data/JSONS",
            feature_type="bonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="mulliken",
            n_jobs=3,
        )

        df = summary_featurize_with_json_bonding.get_df()

        assert isinstance(df, pd.DataFrame)

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
            "edge_COHP",
            "Ionicity_Mull",
        ]

        assert list(df.columns) == expected_cols

    def test_summary_featurize_with_json_antibonding(self):
        summary_featurize_with_json_antibonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="cation-anion",
            path_to_jsons=TestDir / "test_data/Featurizer_test_data/JSONS",
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            charge_type="loewdin",
            n_jobs=3,
        )

        df = summary_featurize_with_json_antibonding.get_df()

        assert isinstance(df, pd.DataFrame)

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
            "edge_COHP",
            "Ionicity_Loew",
        ]

        assert list(df.columns) == expected_cols


class TestBatchCoxxFingerprint:
    def test_fp_cohp_overall(self):
        fp_cohp_overall = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="overall",
            normalize=True,
            tanimoto=True,
            n_jobs=3,
        )
        df = fp_cohp_overall.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(-0.033251, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(-0.013751, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0.046889, abs=1e-05)

    def test_fp_cohp_bonding(self):
        fp_cohp_bonding = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=False,
            tanimoto=True,
            n_jobs=3,
        )
        _ = fp_cohp_bonding.fingerprint_df
        df = fp_cohp_bonding.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0.000017, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0.000000, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0.001532, abs=1e-05)

    def test_fp_cobi(self):
        fp_cobi = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="antibonding",
            normalize=True,
            tanimoto=True,
            fingerprint_for="cobi",
            n_jobs=3,
        )
        _ = fp_cobi.fingerprint_df
        df = fp_cobi.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0, abs=1e-05)

    def test_fp_coop(self):
        fp_coop = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir
            / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=True,
            tanimoto=False,
            fingerprint_for="coop",
            n_jobs=3,
        )
        _ = fp_coop.fingerprint_df
        df = fp_coop.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0, abs=1e-05)


class TestExceptions:
    def test_batch_summary_featurizer_exception(self):
        with pytest.raises(Exception) as err1:  # noqa: PT012, PT011
            self.summary_featurize_with_json_ex = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                / "test_data/Featurizer_test_data/Lobster_calcs_exceptions/1/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json_ex.get_df()

        assert (
            str(err1.value)
            == "COBICAR.lobster or ICOBILIST.lobster file not found in mp-2176"
        )

        with pytest.raises(Exception) as err2:  # noqa: PT012, PT011
            self.summary_featurize_with_json_ex2 = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir
                / "test_data/Featurizer_test_data/Lobster_calcs_exceptions/2/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json_ex2.get_df()

        assert (
            str(err2.value)
            == "COOPCAR.lobster or ICOOPLIST.lobster file not found in mp-1000"
        )

        # COXX exception
        with pytest.raises(Exception) as err3:  # noqa: PT012, PT011
            self.raise_coxx_exception = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/JSONS/"
            )

            _ = self.raise_coxx_exception._featurizecoxx(
                path_to_lobster_calc=self.raise_coxx_exception.path_to_lobster_calcs
            )

        assert (
            str(err3.value)
            == "COHPCAR.lobster or POSCAR or ICOHPLIST.lobster file not found in JSONS"
        )

        # Charges exception
        with pytest.raises(Exception) as err4:  # noqa: PT012, PT011
            self.raise_ch_exception = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/JSONS/"
            )

            _ = self.raise_ch_exception._featurizecharges(
                path_to_lobster_calc=self.raise_ch_exception.path_to_lobster_calcs
            )

        assert str(err4.value) == "CHARGE.lobster or POSCAR not found in JSONS"

        # Fingerprint similarity exception
        with pytest.raises(Exception) as err8:  # noqa: PT012, PT011
            fp_cohp_bonding = BatchCoxxFingerprint(
                path_to_lobster_calcs=TestDir
                / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert (
            str(err8.value)
            == "Cannot compute similarity index. Please set either normalize=True or "
            "tanimoto=True or both to False."
        )
