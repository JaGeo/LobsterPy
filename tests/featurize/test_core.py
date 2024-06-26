from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lobsterpy.featurize.core import (
    DosFingerprint,
    FeaturizeCharges,
    FeaturizeCOXX,
    FeaturizeDoscar,
    FeaturizeLobsterpy,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestFeaturizeLobsterpy:
    def test_featurize_mp1249_json(self):
        featurize_mp1249_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/JSONS/mp-1249.json.gz", bonds="all"
        )
        df = featurize_mp1249_json.get_df(ids="mp-1249")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "mp-1249"

        # Test that all the values in the DataFrame
        assert df.loc["mp-1249", "Icohp_mean_avg"] == pytest.approx(-1.020000, abs=1e-05)
        assert df.loc["mp-1249", "Icohp_mean_max"] == pytest.approx(-1.020000, abs=1e-05)
        assert df.loc["mp-1249", "Icohp_mean_min"] == pytest.approx(-1.020000, abs=1e-05)
        assert df.loc["mp-1249", "Icohp_mean_std"] == pytest.approx(0.000000, abs=1e-05)
        assert df.loc["mp-1249", "Madelung_Mull"] == pytest.approx(-52.000000, abs=1e-05)
        assert df.loc["mp-1249", "bonding_perc_avg"] == pytest.approx(0.978985, abs=1e-05)

    def test_featurize_mp1249_json_ca(self):
        featurize_mp1249_json_ca = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/JSONS/mp-1249.json.gz",
            bonds="cation-anion",
        )
        df = featurize_mp1249_json_ca.get_df(ids="mp-1249")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "mp-1249"

    def test_featurize_mp1958_json(self):
        featurize_mp1958_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/JSONS/mp-1958.json.gz", bonds="all"
        )
        df = featurize_mp1958_json.get_df()

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "mp-1958"

        # Test that all the values in the DataFrame
        assert df.loc["mp-1958", "Icohp_sum_avg"] == pytest.approx(-2.96000, abs=1e-05)
        assert df.loc["mp-1958", "Icohp_sum_max"] == pytest.approx(-2.96000, abs=1e-05)
        assert df.loc["mp-1958", "Icohp_sum_min"] == pytest.approx(-2.96000, abs=1e-05)
        assert df.loc["mp-1958", "Icohp_sum_std"] == pytest.approx(0.000000, abs=1e-05)
        assert df.loc["mp-1958", "Madelung_Loew"] == pytest.approx(-16.68000, abs=1e-05)
        assert df.loc["mp-1958", "antibonding_perc_avg"] == pytest.approx(0.14528, abs=1e-05)

    def test_featurize_mp14652_json(self):
        featurize_mp14652_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/JSONS/mp-14652.json.gz", bonds="all"
        )
        df = featurize_mp14652_json.get_df(ids="mp-14652")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "mp-14652"

        # Test that all the values in the DataFrame
        assert df.loc["mp-14652", "Icohp_mean_std"] == pytest.approx(2.335070, abs=1e-05)
        assert df.loc["mp-14652", "bonding_perc_max"] == pytest.approx(0.889620, abs=1e-05)
        assert df.loc["mp-14652", "bonding_perc_min"] == pytest.approx(0.873420, abs=1e-05)
        assert df.loc["mp-14652", "bonding_perc_std"] == pytest.approx(0.006339, abs=1e-05)
        assert df.loc["mp-14652", "antibonding_perc_min"] == pytest.approx(0.110380, abs=1e-05)
        assert df.loc["mp-14652", "antibonding_perc_max"] == pytest.approx(0.126580, abs=1e-05)
        assert df.loc["mp-14652", "antibonding_perc_std"] == pytest.approx(0.006339, abs=1e-05)

    def test_featurize_mp463(self):
        featurize_mp463 = FeaturizeLobsterpy(
            path_to_lobster_calc=TestDir / "test_data/Featurizer_test_data/Lobster_calcs/mp-463",
            bonds="all",
            orbital_resolved=True,
        )
        df = featurize_mp463.get_df(ids="mp-463")

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
            "Icohp_bndg_orb_mean_avg",
            "Icohp_bndg_orb_mean_max",
            "Icohp_bndg_orb_mean_min",
            "Icohp_bndg_orb_mean_std",
            "Icohp_bndg_orb_sum_avg",
            "Icohp_bndg_orb_sum_max",
            "Icohp_bndg_orb_sum_min",
            "Icohp_bndg_orb_sum_std",
            "bonding_orb_perc_avg",
            "bonding_orb_perc_max",
            "bonding_orb_perc_min",
            "bonding_orb_perc_std",
            "Icohp_antibndg_orb_mean_avg",
            "Icohp_antibndg_orb_mean_max",
            "Icohp_antibndg_orb_mean_min",
            "Icohp_antibndg_orb_mean_std",
            "Icohp_antibndg_orb_sum_avg",
            "Icohp_antibndg_orb_sum_max",
            "Icohp_antibndg_orb_sum_min",
            "Icohp_antibndg_orb_sum_std",
            "antibonding_orb_perc_avg",
            "antibonding_orb_perc_max",
            "antibonding_orb_perc_min",
            "antibonding_orb_perc_std",
            "Madelung_Mull",
            "Madelung_Loew",
        ]

        assert sorted(df.columns) == sorted(expected_cols)

        # Here test now only orbital wise analysis column values
        assert df.loc["mp-463", "Icohp_bndg_orb_mean_avg"] == pytest.approx(-0.181250, abs=1e-05)
        assert df.loc["mp-463", "Icohp_bndg_orb_sum_max"] == pytest.approx(-1.430000, abs=1e-05)
        assert df.loc["mp-463", "bonding_orb_perc_min"] == pytest.approx(0.380000, abs=1e-05)
        assert df.loc["mp-463", "Icohp_antibndg_orb_mean_std"] == pytest.approx(0.059600, abs=1e-05)
        assert df.loc["mp-463", "Icohp_antibndg_orb_sum_avg"] == pytest.approx(-0.714100, abs=1e-05)
        assert df.loc["mp-463", "antibonding_orb_perc_max"] == pytest.approx(0.580000, abs=1e-05)
        assert df.loc["mp-463", "antibonding_orb_perc_std"] == pytest.approx(0.169999, abs=1e-05)

    # Tests for new jsons from atomate2
    def test_featurize_mp66_json(self):
        featurize_mp66_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/Featurizer_test_data/New_JSONS/mp-66.json.gz",
            bonds="all",
        )
        df = featurize_mp66_json.get_df(ids="mp-66")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that all the values in the DataFrame
        assert df.loc["mp-66", "Icohp_mean_avg"] == pytest.approx(-9.59, abs=1e-05)
        assert df.loc["mp-66", "Icohp_sum_max"] == pytest.approx(-38.34, abs=1e-05)
        assert df.loc["mp-66", "Icohp_mean_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "bonding_perc_max"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-66", "bonding_perc_min"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-66", "bonding_perc_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "antibonding_perc_min"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "antibonding_perc_max"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "antibonding_perc_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "Madelung_Mull"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-66", "Madelung_Loew"] == pytest.approx(0, abs=1e-05)

    def test_featurize_mp7000_json(self):
        featurize_mp7000_json = FeaturizeLobsterpy(
            path_to_json=TestDir / "test_data/Featurizer_test_data/New_JSONS/mp-7000.json.gz",
            bonds="cation-anion",
        )
        df = featurize_mp7000_json.get_df(ids="mp-7000")
        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that all the values in the DataFrame
        assert df.loc["mp-7000", "Icohp_mean_avg"] == pytest.approx(-7.98, abs=1e-05)
        assert df.loc["mp-7000", "Icohp_sum_max"] == pytest.approx(-31.90, abs=1e-05)
        assert df.loc["mp-7000", "Icohp_mean_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-7000", "bonding_perc_max"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-7000", "bonding_perc_min"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-7000", "bonding_perc_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-7000", "antibonding_perc_min"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-7000", "antibonding_perc_max"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-7000", "antibonding_perc_std"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-7000", "Madelung_Mull"] == pytest.approx(-163.37, abs=1e-05)
        assert df.loc["mp-7000", "Madelung_Loew"] == pytest.approx(-99.01, abs=1e-05)

    def test_featurize_csh_madelung(self):
        featurize_csh_madelung = FeaturizeLobsterpy(path_to_lobster_calc=TestDir / "test_data/CsH/", bonds="all")
        df = featurize_csh_madelung.get_df()

        assert np.isnan(df.loc["CsH", "Madelung_Mull"])
        assert np.isnan(df.loc["CsH", "Madelung_Loew"])


class TestFeaturizeCOXX:
    def test_featurize_nacl_coxx(self):
        featurize_nacl_coxx = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
            path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
            path_to_structure=TestDir / "test_data/NaCl/POSCAR.gz",
            feature_type="overall",
            e_range=[-5, 0],
        )
        df = featurize_nacl_coxx.get_summarized_coxx_df(ids="NaCl")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "NaCl"

        # Test that all the values in the DataFrame
        assert df.loc["NaCl", "bnd_wICOHP"] == pytest.approx(97.233878, abs=1e-05)
        assert df.loc["NaCl", "antibnd_wICOHP"] == pytest.approx(2.766122, abs=1e-05)
        assert df.loc["NaCl", "w_ICOHP"] == pytest.approx(-0.150558, abs=1e-05)

        assert df.loc["NaCl", "EIN_ICOHP"] == pytest.approx(27.843536, abs=1e-05)
        assert df.loc["NaCl", "center_COHP"] == pytest.approx(-4.96241, abs=1e-05)
        assert df.loc["NaCl", "width_COHP"] == pytest.approx(8.881784e-16, abs=1e-05)
        assert df.loc["NaCl", "skewness_COHP"] == pytest.approx(1, abs=1e-05)
        assert df.loc["NaCl", "kurtosis_COHP"] == pytest.approx(1, abs=1e-05)

        # test summary features using label list
        df1 = featurize_nacl_coxx.get_summarized_coxx_df(label_list=["2", "3"])
        assert df.loc["NaCl", "center_COHP"] != df1.loc["NaCl", "center_COHP"]

        df_fp = featurize_nacl_coxx.get_coxx_fingerprint_df(n_bins=20000)

        fingerprint = df_fp.loc["NaCl", "COXX_FP"]

        assert fingerprint.n_bins != 20000

        df_fp1 = featurize_nacl_coxx.get_coxx_fingerprint_df(binning=False)

        fingerprint = df_fp1.loc["NaCl", "COXX_FP"]

        assert fingerprint.n_bins == 401

        df_fp2 = featurize_nacl_coxx.get_coxx_fingerprint_df(label_list=["3", "5"])

        fingerprint_label = df_fp2.loc["NaCl", "COXX_FP"]

        assert fingerprint.__str__() != fingerprint_label.__str__()

        # test moment features using label and orbital list
        label_list = ["21", "23", "24", "27", "28", "30"]
        (coxx_c, coxx_w, coxx_s, coxx_k, coxx_e) = FeaturizeCOXX._calc_moment_features(
            complete_coxx_obj=featurize_nacl_coxx.completecoxx,
            e_range=[None, 0],
            feature_type=featurize_nacl_coxx.feature_type,
            label_list=label_list,
            orbital="3px-2px",
        )

        assert coxx_c == pytest.approx(-5.8369822833428735, abs=1e-05)
        assert coxx_w == pytest.approx(0.4523359129845927, abs=1e-05)
        assert coxx_s == pytest.approx(-0.15433775869720565, abs=1e-05)
        assert coxx_k == pytest.approx(1.6916251153182233, abs=1e-05)
        assert coxx_e == pytest.approx(-5.26316, abs=1e-05)

    def test_featurize_cdf_coxx(self):
        featurize_cdf_coxx = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "test_data/CdF/COHPCAR.lobster.gz",
            path_to_icoxxlist=TestDir / "test_data/CdF/ICOHPLIST.lobster.gz",
            path_to_structure=TestDir / "test_data/CdF/POSCAR.gz",
            feature_type="bonding",
            e_range=[-5, 0],
        )
        df = featurize_cdf_coxx.get_summarized_coxx_df()

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "CdF"

        # Test that all the values in the DataFrame
        assert df.loc["CdF", "bnd_wICOHP"] == pytest.approx(81.112732, abs=1e-05)
        assert df.loc["CdF", "antibnd_wICOHP"] == pytest.approx(18.887268, abs=1e-05)
        assert df.loc["CdF", "w_ICOHP"] == pytest.approx(-0.198235, abs=1e-05)

        assert df.loc["CdF", "EIN_ICOHP"] == pytest.approx(18.76634, abs=1e-05)
        assert df.loc["CdF", "center_COHP"] == pytest.approx(-4.748383, abs=1e-05)
        assert df.loc["CdF", "width_COHP"] == pytest.approx(0.157761, abs=1e-05)
        assert df.loc["CdF", "skewness_COHP"] == pytest.approx(0.910094, abs=1e-05)
        assert df.loc["CdF", "kurtosis_COHP"] == pytest.approx(2.866611, abs=1e-05)
        assert df.loc["CdF", "edge_COHP"] == pytest.approx(-4.96241, abs=1e-05)

        # test using label list
        df1 = featurize_cdf_coxx.get_summarized_coxx_df(label_list=["2", "3", "30"])
        assert df.loc["CdF", "center_COHP"] != df1.loc["CdF", "center_COHP"]

        df_fp = featurize_cdf_coxx.get_coxx_fingerprint_df(n_bins=20000)

        fingerprint = df_fp.loc["CdF", "COXX_FP"]

        assert fingerprint.n_bins != 20000

        df_fp1 = featurize_cdf_coxx.get_coxx_fingerprint_df(binning=False)

        fingerprint = df_fp1.loc["CdF", "COXX_FP"]

        assert fingerprint.n_bins == 401

        df_fp2 = featurize_cdf_coxx.get_coxx_fingerprint_df(label_list=["3", "5"])

        fingerprint_label = df_fp2.loc["CdF", "COXX_FP"]

        assert fingerprint.__str__() != fingerprint_label.__str__()

        # test moment features using label and orbital list
        label_list = ["25", "32", "35", "36", "57", "58", "61", "68"]
        (coxx_c, coxx_w, coxx_s, coxx_k, coxx_e) = FeaturizeCOXX._calc_moment_features(
            complete_coxx_obj=featurize_cdf_coxx.completecoxx,
            e_range=featurize_cdf_coxx.e_range,
            feature_type=featurize_cdf_coxx.feature_type,
            label_list=label_list,
            orbital="2py-5s",
        )

        assert coxx_c == pytest.approx(-4.609186587118677, abs=1e-05)
        assert coxx_w == pytest.approx(0.2416386973020571, abs=1e-05)
        assert coxx_s == pytest.approx(0.7943687717963166, abs=1e-05)
        assert coxx_k == pytest.approx(3.196850488714957, abs=1e-05)
        assert coxx_e == pytest.approx(-4.96241, abs=1e-05)

    def test_featurize_k3sb_coxx(self):
        featurize_k3sb_coxx = FeaturizeCOXX(
            path_to_coxxcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
            path_to_icoxxlist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
            path_to_structure=TestDir / "test_data/K3Sb/POSCAR.gz",
            feature_type="antibonding",
            e_range=[None, None],
        )
        df = featurize_k3sb_coxx.get_summarized_coxx_df(ids="K3Sb")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

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
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "K3Sb"

        # Test that all the values in the DataFrame
        assert df.loc["K3Sb", "bnd_wICOHP"] == pytest.approx(97.019044, abs=1e-05)
        assert df.loc["K3Sb", "antibnd_wICOHP"] == pytest.approx(2.980956, abs=1e-05)
        assert df.loc["K3Sb", "w_ICOHP"] == pytest.approx(-0.318218, abs=1e-05)

        assert df.loc["K3Sb", "EIN_ICOHP"] == pytest.approx(11.597595, abs=1e-05)
        assert df.loc["K3Sb", "center_COHP"] == pytest.approx(2.786451, abs=1e-05)
        assert df.loc["K3Sb", "width_COHP"] == pytest.approx(0.765492, abs=1e-05)
        assert df.loc["K3Sb", "skewness_COHP"] == pytest.approx(-1.52563, abs=1e-05)
        assert df.loc["K3Sb", "kurtosis_COHP"] == pytest.approx(6.829327, abs=1e-05)
        assert df.loc["K3Sb", "edge_COHP"] == pytest.approx(3.1916, abs=1e-05)

        # test moment features using label and orbital list
        label_list = [
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "21",
            "22",
            "23",
            "24",
            "25",
            "26",
            "27",
            "28",
        ]
        (coxx_c, coxx_w, coxx_s, coxx_k, coxx_e) = FeaturizeCOXX._calc_moment_features(
            complete_coxx_obj=featurize_k3sb_coxx.completecoxx,
            e_range=featurize_k3sb_coxx.e_range,
            feature_type=featurize_k3sb_coxx.feature_type,
            label_list=label_list,
            orbital="5s-4s",
        )

        assert coxx_c == pytest.approx(1.4416312688111554, abs=1e-05)
        assert coxx_w == pytest.approx(3.4768865885682962, abs=1e-05)
        assert coxx_s == pytest.approx(-7.434008678770313, abs=1e-05)
        assert coxx_k == pytest.approx(70.35977484462431, abs=1e-05)
        assert coxx_e == pytest.approx(-31.56499, abs=1e-05)


class TestFeaturizeCharges:
    def test_featurize_c_charge(self):
        featurize_c_charge = FeaturizeCharges(
            path_to_structure=TestDir / "test_data/C/POSCAR.gz",
            path_to_charge=TestDir / "test_data/C/CHARGE.lobster.gz",
            charge_type="mulliken",
        )
        df = featurize_c_charge.get_df(ids="C")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Mull",
        ]
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "C"

        # Test that all the values in the DataFrame
        assert df.loc["C", "Ionicity_Mull"] == pytest.approx(0.0, abs=1e-05)

    def test_featurize_cdf_charge(self):
        featurize_cdf_charge = FeaturizeCharges(
            path_to_structure=TestDir / "test_data/CdF/POSCAR.gz",
            path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster.gz",
            charge_type="mulliken",
        )
        df = featurize_cdf_charge.get_df(ids="CdF")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Mull",
        ]
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "CdF"

        # Test that all the values in the DataFrame
        assert df.loc["CdF", "Ionicity_Mull"] == pytest.approx(0.788333, abs=1e-05)

    def test_featurize_k3sb_charge(self):
        featurize_k3sb_charge = FeaturizeCharges(
            path_to_structure=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            charge_type="loewdin",
        )
        df = featurize_k3sb_charge.get_df(ids="K3Sb")

        # Test that the function returns a pandas DataFrame
        assert isinstance(df, pd.DataFrame)

        # Test that the DataFrame has the expected columns
        expected_cols = [
            "Ionicity_Loew",
        ]
        assert sorted(df.columns) == sorted(expected_cols)

        # Test that the DataFrame has the expected index
        assert df.index[0] == "K3Sb"

        # Test that all the values in the DataFrame
        assert df.loc["K3Sb", "Ionicity_Loew"] == pytest.approx(0.563333, abs=1e-05)


class TestExceptions:
    def test_lobsterpy_featurize_exception(self):
        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_mp1249_json = FeaturizeLobsterpy(path_to_json=None, path_to_lobster_calc=None, bonds="all")

            _ = self.featurize_mp1249_json.get_df()

        assert str(err.value) == "Please provide either path to lightweight lobster jsons or path to lobster calc"

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_mp1249_json = FeaturizeLobsterpy(
                path_to_json=None, path_to_lobster_calc=TestDir, bonds="all"
            )

            _ = self.featurize_mp1249_json.get_df()

        assert (
            str(err.value) == "Files ['POSCAR', 'COHPCAR.lobster', 'ICOHPLIST.lobster', 'CHARGE.lobster'] "
            "not found in ..."
        )

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_csh_cation_anion = FeaturizeLobsterpy(
                path_to_lobster_calc=TestDir / "test_data/CsH/", bonds="cation-anion"
            )

            _ = self.featurize_csh_cation_anion.get_df()

        assert str(err.value) == "No cation-anion bonds detected for CsH structure. Please switch to `all` bonds mode"

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_c_cation_anion = FeaturizeLobsterpy(
                path_to_lobster_calc=TestDir / "test_data/C/", bonds="cation-anion"
            )

            _ = self.featurize_c_cation_anion.get_df()

        assert str(err.value) == "No cation-anion bonds detected for C structure. Please switch to `all` bonds mode"

    def test_featurize_charges(self):
        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_cdf_charge = FeaturizeCharges(
                path_to_structure=TestDir / "test_data/CdF/POSCAR.gz",
                path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster.gz",
                charge_type="Mull",
            )

            _ = self.featurize_cdf_charge.get_df()

        assert str(err.value) == "Please check the requested charge_type. Possible options are `mulliken` or `loewdin`"

    def test_featurize_coxx(self):
        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_coxx = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
                path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
                path_to_structure=TestDir / "test_data/NaCl/POSCAR.gz",
                feature_type="summed",
                e_range=[None, None],
            )

            _ = self.featurize_coxx.get_summarized_coxx_df(ids="exception")

        assert (
            str(err.value)
            == "Please recheck feature type requested argument. Possible options are bonding/antibonding/overall"
        )

        with pytest.raises(Exception) as err2:  # noqa: PT012, PT011
            self.featurize_coxx = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
                path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
                path_to_structure=TestDir / "test_data/NaCl/POSCAR.gz",
                feature_type="bonding",
                e_range=[None, None],
            )

            _ = self.featurize_coxx.get_coxx_fingerprint_df(spin_type="-1")

        assert str(err2.value) == "Check the spin_type argument.Possible options are summed/up/down"

        with pytest.raises(Exception) as err3:  # noqa: PT012, PT011
            self.featurize_coxx = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "test_data/NaSi/COHPCAR.lobster.gz",
                path_to_icoxxlist=TestDir / "test_data/NaSi/ICOHPLIST.lobster.gz",
                path_to_structure=TestDir / "test_data/NaSi/POSCAR.gz",
                feature_type="bonding",
                e_range=[-5, 0],
                are_cobis=True,
                are_coops=True,
            )

            _ = self.featurize_coxx.get_coxx_fingerprint_df()

        assert str(err3.value) == "You cannot have info about COOPs and COBIs in the same file."

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_nacl_coxx = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
                path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
                path_to_structure=TestDir / "test_data/NaCl/POSCAR.gz",
                feature_type="antibond",
                e_range=[-5, 0],
            )

            _ = self.featurize_nacl_coxx.get_coxx_fingerprint_df(ids="NACL")

        assert (
            str(err.value) == "Please recheck fingerprint type requested argument. "
            "Possible options are bonding/antibonding/overall"
        )

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            self.featurize_nacl_coxx = FeaturizeCOXX(
                path_to_coxxcar=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
                path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
                path_to_structure=TestDir / "test_data/NaCl/POSCAR.gz",
                feature_type="antibonding",
                e_range=[-5, 0],
            )

            _ = self.featurize_nacl_coxx.get_coxx_fingerprint_df(spin_type="down")

        assert str(err.value) == "LOBSTER calculation is non-spin polarized. Please switch spin_type to `up`"


class TestFeaturizeDoscar:
    def test_featurize_nacl_dos(self):
        feat_dos = FeaturizeDoscar(
            path_to_doscar=TestDir / "test_data/NaCl_comp_range/DOSCAR.LSO.lobster.gz",
            path_to_structure=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            e_range=[-10, 0],
        )

        assert feat_dos.get_df().index[0] == "NaCl_comp_range"

        assert feat_dos.get_df(ids="NaCl").index[0] == "NaCl"

        df = feat_dos.get_df(ids="NaCl")

        # Test that all the values in the DataFrame
        assert df.loc["NaCl", "s_band_center"] == pytest.approx(-1.3175, abs=1e-05)
        assert df.loc["NaCl", "s_band_width"] == pytest.approx(0.3207, abs=1e-05)
        assert df.loc["NaCl", "s_band_skew"] == pytest.approx(1.1887, abs=1e-05)
        assert df.loc["NaCl", "s_band_kurtosis"] == pytest.approx(4.0563, abs=1e-05)
        assert df.loc["NaCl", "s_band_upperband_edge"] == pytest.approx(-1.3941, abs=1e-05)

        assert df.loc["NaCl", "p_band_center"] == pytest.approx(-0.7484, abs=1e-05)
        assert df.loc["NaCl", "p_band_width"] == pytest.approx(0.4815, abs=1e-05)
        assert df.loc["NaCl", "p_band_skew"] == pytest.approx(-0.5047, abs=1e-05)
        assert df.loc["NaCl", "p_band_kurtosis"] == pytest.approx(1.8789, abs=1e-05)
        assert df.loc["NaCl", "p_band_upperband_edge"] == pytest.approx(-0.2343, abs=1e-05)

    def test_featurize_k3sb_dos(self):
        feat_dos = FeaturizeDoscar(
            path_to_doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
            path_to_structure=TestDir / "test_data/K3Sb/POSCAR.gz",
            add_element_dos_moments=True,
            e_range=None,
        )

        assert feat_dos.get_df().index[0] == "K3Sb"

        df = feat_dos.get_df(ids="K3Sb")

        # Test that all the orbital moment values in the DataFrame
        assert df.loc["K3Sb", "s_band_center"] == pytest.approx(-13.3722, abs=1e-05)
        assert df.loc["K3Sb", "s_band_width"] == pytest.approx(15.5141, abs=1e-05)
        assert df.loc["K3Sb", "s_band_skew"] == pytest.approx(-0.1718, abs=1e-05)
        assert df.loc["K3Sb", "s_band_kurtosis"] == pytest.approx(1.1523, abs=1e-05)
        assert df.loc["K3Sb", "s_band_upperband_edge"] == pytest.approx(-31.5650, abs=1e-05)

        assert df.loc["K3Sb", "p_band_center"] == pytest.approx(-10.7245, abs=1e-05)
        assert df.loc["K3Sb", "p_band_width"] == pytest.approx(6.4334, abs=1e-05)
        assert df.loc["K3Sb", "p_band_skew"] == pytest.approx(1.1553, abs=1e-05)
        assert df.loc["K3Sb", "p_band_kurtosis"] == pytest.approx(2.5024, abs=1e-05)
        assert df.loc["K3Sb", "p_band_upperband_edge"] == pytest.approx(-14.0357, abs=1e-05)

        # Test Element orbital moment features

        assert df.loc["K3Sb", "K_s_band_center"] == pytest.approx(-14.6044, abs=1e-05)
        assert df.loc["K3Sb", "K_s_band_width"] == pytest.approx(16.3903, abs=1e-05)
        assert df.loc["K3Sb", "K_s_band_skew"] == pytest.approx(0.0309, abs=1e-05)
        assert df.loc["K3Sb", "K_s_band_kurtosis"] == pytest.approx(1.0352, abs=1e-05)
        assert df.loc["K3Sb", "K_s_band_upperband_edge"] == pytest.approx(-31.5650, abs=1e-05)

        assert df.loc["K3Sb", "Sb_s_band_center"] == pytest.approx(-5.9795, abs=1e-05)
        assert df.loc["K3Sb", "Sb_s_band_width"] == pytest.approx(3.0357, abs=1e-05)
        assert df.loc["K3Sb", "Sb_s_band_skew"] == pytest.approx(1.9959, abs=1e-05)
        assert df.loc["K3Sb", "Sb_s_band_kurtosis"] == pytest.approx(8.4750, abs=1e-05)
        assert df.loc["K3Sb", "Sb_s_band_upperband_edge"] == pytest.approx(-6.9954, abs=1e-05)

        assert df.loc["K3Sb", "K_p_band_center"] == pytest.approx(-14.3425, abs=1e-05)
        assert df.loc["K3Sb", "K_p_band_width"] == pytest.approx(1.1348, abs=1e-05)
        assert df.loc["K3Sb", "K_p_band_skew"] == pytest.approx(10.7021, abs=1e-05)
        assert df.loc["K3Sb", "K_p_band_kurtosis"] == pytest.approx(147.5038, abs=1e-05)
        assert df.loc["K3Sb", "K_p_band_upperband_edge"] == pytest.approx(-14.0357, abs=1e-05)

        assert df.loc["K3Sb", "Sb_p_band_center"] == pytest.approx(0.1300, abs=1e-05)
        assert df.loc["K3Sb", "Sb_p_band_width"] == pytest.approx(2.1458, abs=1e-05)
        assert df.loc["K3Sb", "Sb_p_band_skew"] == pytest.approx(-6.1882, abs=1e-05)
        assert df.loc["K3Sb", "Sb_p_band_kurtosis"] == pytest.approx(84.5825, abs=1e-05)
        assert df.loc["K3Sb", "Sb_p_band_upperband_edge"] == pytest.approx(-0.0775, abs=1e-05)

        # Test for the case where e_range is set to None and trying to get fingerprint
        df_fp = feat_dos.get_fingerprint_df()
        assert isinstance(df_fp, pd.DataFrame)
        assert isinstance(df_fp.loc["K3Sb", "DOS_FP"], DosFingerprint)
