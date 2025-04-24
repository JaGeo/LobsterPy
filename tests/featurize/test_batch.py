import re
from pathlib import Path

import pandas as pd
import pytest
from numpy import count_nonzero, unique
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.electronic_structure.dos import DosFingerprint

from lobsterpy.featurize.batch import (
    BatchCoxxFingerprint,
    BatchDosFeaturizer,
    BatchIcoxxlistFeaturizer,
    BatchStructureGraphs,
    BatchSummaryFeaturizer,
)
from lobsterpy.featurize.core import CoxxFingerprint

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestBatchSummaryFeaturizer:
    def test_summary_featurize_with_json(self):
        summary_featurize_with_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert sorted(df.columns) == sorted(expected_cols)

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert sorted(df.index) == sorted(expected_index)

    def test_summary_featurize_with_no_bonds(self):
        summary_featurize = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/No_bonds_cases",
            bonds="all",
            path_to_jsons=None,
            feature_type="antibonding",
            include_cobi_data=False,
            include_coop_data=False,
            e_range=[-15, 0],
            n_jobs=2,
        )

        df = summary_featurize.get_df()

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

        assert sorted(df.columns) == sorted(expected_cols)

        expected_index = ["mp-111", "mp-23155"]

        assert sorted(df.index) == sorted(expected_index)

    def test_summary_featurize_orbitalwise(self):
        summary_featurize_without_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            bonds="all",
            include_cobi_data=False,
            include_coop_data=False,
            orbital_resolved=True,
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

        assert sorted(df.columns) == sorted(expected_cols)

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert sorted(df.index) == sorted(expected_index)

    def test_summary_featurize_without_json(self):
        summary_featurize_without_json = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert sorted(df.columns) == sorted(expected_cols)

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert sorted(df.index) == sorted(expected_index)

    def test_summary_featurize_with_json_overall(self):
        summary_featurize_with_json_overall = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert sorted(df.columns) == sorted(expected_cols)

        expected_index = ["mp-1000", "mp-2176", "mp-463"]

        assert sorted(df.index) == sorted(expected_index)

    def test_summary_featurize_with_json_bonding(self):
        summary_featurize_with_json_bonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert sorted(df.columns) == sorted(expected_cols)

    def test_summary_featurize_with_json_antibonding(self):
        summary_featurize_with_json_antibonding = BatchSummaryFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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

        assert sorted(df.columns) == sorted(expected_cols)


class TestBatchCoxxFingerprint:
    def test_fp_cohp_overall(self):
        fp_cohp_overall = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=False,
            tanimoto=True,
            n_jobs=3,
        )
        df_fp_cohp = fp_cohp_bonding.fingerprint_df
        for fp in df_fp_cohp["COXX_FP"]:
            assert isinstance(fp, CoxxFingerprint)
            assert fp.fp_type == "bonding"
        df = fp_cohp_bonding.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0.000017, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0.000000, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0.001532, abs=1e-05)

    def test_fp_cobi(self):
        fp_cobi = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="antibonding",
            normalize=True,
            tanimoto=True,
            fingerprint_for="cobi",
            n_jobs=3,
        )
        df_fp_cobi = fp_cobi.fingerprint_df
        for fp in df_fp_cobi["COXX_FP"]:
            assert isinstance(fp, CoxxFingerprint)
            assert fp.fp_type == "antibonding"

        df = fp_cobi.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0, abs=1e-05)

    def test_fp_coop(self):
        fp_coop = BatchCoxxFingerprint(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            e_range=[-15, 0],
            feature_type="bonding",
            normalize=True,
            tanimoto=False,
            fingerprint_for="coop",
            n_jobs=3,
        )
        df_fp_coop = fp_coop.fingerprint_df
        for fp in df_fp_coop["COXX_FP"]:
            assert isinstance(fp, CoxxFingerprint)
            assert fp.fp_type == "bonding"

        df = fp_coop.get_similarity_matrix_df()

        assert df.loc["mp-463", "mp-1000"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-2176"] == pytest.approx(0, abs=1e-05)
        assert df.loc["mp-463", "mp-463"] == pytest.approx(1, abs=1e-05)
        assert df.loc["mp-1000", "mp-2176"] == pytest.approx(0, abs=1e-05)


class TestBatchStructureGraphs:
    def test_batch_structure_graphs_all_bonds(self):
        batch_sg = BatchStructureGraphs(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            which_bonds="all",
            n_jobs=3,
        )

        df = batch_sg.get_df()

        assert isinstance(df, pd.DataFrame)
        for graph_obj in df["structure_graph"]:
            assert isinstance(graph_obj, StructureGraph)
        assert len(df.index) == 3

    def test_batch_structure_graphs_cation_anion_bonds(self):
        batch_sg = BatchStructureGraphs(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            which_bonds="cation-anion",
            n_jobs=3,
        )

        df = batch_sg.get_df()

        assert isinstance(df, pd.DataFrame)
        for graph_obj in df["structure_graph"]:
            assert isinstance(graph_obj, StructureGraph)
        assert len(df.index) == 3


class TestBatchDosFeaturizer:
    def test_batch_dos_featurizer_non_lso(self):
        batch_dos = BatchDosFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            use_lso_dos=False,
            e_range=[-5, 0],
            fingerprint_type="p",
            n_bins=100,
            n_jobs=3,
        )

        df_moments = batch_dos.get_df()
        assert isinstance(df_moments, pd.DataFrame)

        expected_cols = [
            "s_band_center",
            "s_band_width",
            "s_band_skew",
            "s_band_kurtosis",
            "s_band_upperband_edge",
            "p_band_center",
            "p_band_width",
            "p_band_skew",
            "p_band_kurtosis",
            "p_band_upperband_edge",
            "d_band_center",
            "d_band_width",
            "d_band_skew",
            "d_band_kurtosis",
            "d_band_upperband_edge",
        ]

        assert sorted(df_moments.columns) == sorted(expected_cols)

        df_fp = batch_dos.get_fingerprints_df()
        assert isinstance(df_fp, pd.DataFrame)

        for dos_fp in df_fp["DOS_FP"]:
            assert isinstance(dos_fp, DosFingerprint)
            assert dos_fp.fp_type == "p"
            assert dos_fp.n_bins == 100

    def test_batch_dos_featurizer_lso(self):
        batch_dos_lso = BatchDosFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            add_element_dos_moments=True,
            use_lso_dos=True,
            e_range=[-5, 0],
            fingerprint_type="summed_pdos",
            n_bins=256,
            n_jobs=3,
        )

        df_moments_lso = batch_dos_lso.get_df()
        assert isinstance(df_moments_lso, pd.DataFrame)

        expected_cols = [
            "s_band_center",
            "s_band_width",
            "s_band_skew",
            "s_band_kurtosis",
            "s_band_upperband_edge",
            "Ba_s_band_center",
            "Ba_s_band_width",
            "Ba_s_band_skew",
            "Ba_s_band_kurtosis",
            "Ba_s_band_upperband_edge",
            "Te_s_band_center",
            "Te_s_band_width",
            "Te_s_band_skew",
            "Te_s_band_kurtosis",
            "Te_s_band_upperband_edge",
            "p_band_center",
            "p_band_width",
            "p_band_skew",
            "p_band_kurtosis",
            "p_band_upperband_edge",
            "Ba_p_band_center",
            "Ba_p_band_width",
            "Ba_p_band_skew",
            "Ba_p_band_kurtosis",
            "Ba_p_band_upperband_edge",
            "Te_p_band_center",
            "Te_p_band_width",
            "Te_p_band_skew",
            "Te_p_band_kurtosis",
            "Te_p_band_upperband_edge",
            "Zn_s_band_center",
            "Zn_s_band_width",
            "Zn_s_band_skew",
            "Zn_s_band_kurtosis",
            "Zn_s_band_upperband_edge",
            "d_band_center",
            "d_band_width",
            "d_band_skew",
            "d_band_kurtosis",
            "d_band_upperband_edge",
            "Zn_d_band_center",
            "Zn_d_band_width",
            "Zn_d_band_skew",
            "Zn_d_band_kurtosis",
            "Zn_d_band_upperband_edge",
            "K_s_band_center",
            "K_s_band_width",
            "K_s_band_skew",
            "K_s_band_kurtosis",
            "K_s_band_upperband_edge",
            "F_s_band_center",
            "F_s_band_width",
            "F_s_band_skew",
            "F_s_band_kurtosis",
            "F_s_band_upperband_edge",
            "K_p_band_center",
            "K_p_band_width",
            "K_p_band_skew",
            "K_p_band_kurtosis",
            "K_p_band_upperband_edge",
            "F_p_band_center",
            "F_p_band_width",
            "F_p_band_skew",
            "F_p_band_kurtosis",
            "F_p_band_upperband_edge",
        ]

        assert sorted(df_moments_lso.columns) == sorted(expected_cols)

        df_fp_lso = batch_dos_lso.get_fingerprints_df()
        assert isinstance(df_fp_lso, pd.DataFrame)

        for dos_fp in df_fp_lso["DOS_FP"]:
            assert isinstance(dos_fp, DosFingerprint)
            assert dos_fp.fp_type == "summed_pdos"
            assert dos_fp.n_bins == 256


class TestBatchIcoxxlistFeaturizer:
    @pytest.mark.parametrize(
        ("normalization", "bwdf_df_type", "sorted_dists_mode"),
        [
            ("formula_units", "binned", "negative"),
            ("area", "stats", "negative"),
            ("counts", "stats", "negative"),
            ("none", "binned", "negative"),
            ("formula_units", "sorted_dists", "negative"),
            ("formula_units", "sorted_bwdf", "negative"),
        ],
    )
    def test_batch_icohplist_featurizer(self, normalization, bwdf_df_type, sorted_dists_mode):
        batch_icohp = BatchIcoxxlistFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            n_jobs=3,
            bin_width=0.5,
            normalization=normalization,
            bwdf_df_type=bwdf_df_type,
            sorted_dists_mode=sorted_dists_mode,
        )

        df_icohp = batch_icohp.get_df()
        expected_index = ["mp-1000", "mp-2176", "mp-463"]
        assert isinstance(df_icohp, pd.DataFrame)
        assert sorted(df_icohp.index) == sorted(expected_index)
        if bwdf_df_type != "sorted_dists":
            assert not (df_icohp >= 0).all().all()  # Test if all values are < zero (icohps are read correctly)
        if bwdf_df_type == "binned":
            expected_cols = [
                "bwdf_0.0-0.5",
                "bwdf_0.5-1.0",
                "bwdf_1.0-1.5",
                "bwdf_1.5-2.0",
                "bwdf_2.0-2.5",
                "bwdf_2.5-3.0",
                "bwdf_3.0-3.5",
                "bwdf_3.5-4.0",
                "bwdf_4.0-4.5",
                "bwdf_4.5-5.0",
                "bwdf_5.0-5.5",
                "bwdf_5.5-6.0",
            ]
            assert sorted(df_icohp.columns) == sorted(expected_cols)

        elif bwdf_df_type == "stats":
            stats_df_expected_columns = [
                "bwdf_sum",
                "bwdf_mean",
                "bwdf_std",
                "bwdf_min",
                "bwdf_max",
                "bwdf_skew",
                "bwdf_kurtosis",
                "bwdf_w_mean",
                "bwdf_w_std",
            ]
            assert sorted(df_icohp.columns) == sorted(stats_df_expected_columns)

        elif bwdf_df_type == "sorted_bwdf":
            assert not df_icohp.isna().to_numpy().any()
            unique_values, counts = unique(df_icohp.values, return_counts=True)
            assert dict(zip(unique_values, counts))[0.0] == 3
            assert set(df_icohp.columns) == {f"bwdf_at_dist{c_idx}" for c_idx in range(4)}

        else:  # sorted_dists
            assert not df_icohp.isna().to_numpy().any()
            unique_values, counts = unique(df_icohp.values, return_counts=True)
            assert dict(zip(unique_values, counts))[0.0] == 3
            assert set(df_icohp.columns) == {f"dist_at_{sorted_dists_mode[:3]}_bwdf{c_idx}" for c_idx in range(4)}

    @pytest.mark.parametrize(
        ("normalization", "bwdf_df_type", "sorted_dists_mode"),
        [
            ("formula_units", "stats", "negative"),
            ("area", "binned", "negative"),
            ("counts", "binned", "negative"),
            ("none", "stats", "negative"),
            ("area", "sorted_dists", "positive"),
            ("area", "sorted_dists", "negative"),
            ("area", "sorted_bwdf", "negative"),
        ],
    )
    def test_batch_icobilist_featurizer(self, normalization, bwdf_df_type, sorted_dists_mode):
        batch_icobi = BatchIcoxxlistFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            n_jobs=3,
            bin_width=0.5,
            read_icobis=True,
            normalization=normalization,
            bwdf_df_type=bwdf_df_type,
            sorted_dists_mode=sorted_dists_mode,
        )

        df_icobi = batch_icobi.get_df()
        expected_index = ["mp-1000", "mp-2176", "mp-463"]
        # Test if all values are above zero > icobis are read
        assert (df_icobi >= 0).all().all()  # check if all values are above zero
        assert isinstance(df_icobi, pd.DataFrame)
        assert sorted(df_icobi.index) == sorted(expected_index)

        if bwdf_df_type == "binned":
            expected_cols = [
                "bwdf_0.0-0.5",
                "bwdf_0.5-1.0",
                "bwdf_1.0-1.5",
                "bwdf_1.5-2.0",
                "bwdf_2.0-2.5",
                "bwdf_2.5-3.0",
                "bwdf_3.0-3.5",
                "bwdf_3.5-4.0",
                "bwdf_4.0-4.5",
                "bwdf_4.5-5.0",
                "bwdf_5.0-5.5",
                "bwdf_5.5-6.0",
            ]
            assert sorted(df_icobi.columns) == sorted(expected_cols)

        elif bwdf_df_type == "stats":
            stats_df_expected_columns = [
                "bwdf_sum",
                "bwdf_mean",
                "bwdf_std",
                "bwdf_min",
                "bwdf_max",
                "bwdf_skew",
                "bwdf_kurtosis",
                "bwdf_w_mean",
                "bwdf_w_std",
            ]
            assert sorted(df_icobi.columns) == sorted(stats_df_expected_columns)

        elif bwdf_df_type == "sorted_bwdf":
            assert not df_icobi.isna().to_numpy().any()
            unique_values = unique(df_icobi.values)
            assert 0.0 not in unique_values
            assert set(df_icobi.columns) == {f"bwdf_at_dist{c_idx}" for c_idx in range(2)}

        else:  # sorted_dists
            assert not df_icobi.isna().to_numpy().any()
            unique_values = unique(df_icobi.values)
            if sorted_dists_mode == "positive":
                assert 0.0 not in unique_values
                assert set(df_icobi.columns) == {f"dist_at_{sorted_dists_mode[:3]}_bwdf{c_idx}" for c_idx in range(2)}
            else:
                assert df_icobi.empty

    @pytest.mark.parametrize(
        ("normalization", "bwdf_df_type", "sorted_dists_mode"),
        [
            ("formula_units", "binned", "positive"),
            ("area", "stats", "positive"),
            ("counts", "stats", "positive"),
            ("none", "binned", "positive"),
            ("none", "sorted_dists", "positive"),
            ("none", "sorted_dists", "negative"),
            ("none", "sorted_bwdf", "positive"),
        ],
    )
    def test_batch_icooplist_featurizer(self, normalization, bwdf_df_type, sorted_dists_mode):
        batch_icoop = BatchIcoxxlistFeaturizer(
            path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
            n_jobs=3,
            bin_width=0.5,
            read_icoops=True,
            normalization=normalization,
            bwdf_df_type=bwdf_df_type,
            sorted_dists_mode=sorted_dists_mode,
        )

        df_icoop = batch_icoop.get_df()
        expected_index = ["mp-1000", "mp-2176", "mp-463"]
        if bwdf_df_type != "sorted_dists":
            assert not (df_icoop >= 0).all().all()  # Test if all values are above zero > icobis are read
        assert isinstance(df_icoop, pd.DataFrame)
        assert sorted(df_icoop.index) == sorted(expected_index)

        if bwdf_df_type == "binned":
            expected_cols = [
                "bwdf_0.0-0.5",
                "bwdf_0.5-1.0",
                "bwdf_1.0-1.5",
                "bwdf_1.5-2.0",
                "bwdf_2.0-2.5",
                "bwdf_2.5-3.0",
                "bwdf_3.0-3.5",
                "bwdf_3.5-4.0",
                "bwdf_4.0-4.5",
                "bwdf_4.5-5.0",
                "bwdf_5.0-5.5",
                "bwdf_5.5-6.0",
            ]
            assert sorted(df_icoop.columns) == sorted(expected_cols)

        elif bwdf_df_type == "stats":
            stats_df_expected_columns = [
                "bwdf_sum",
                "bwdf_mean",
                "bwdf_std",
                "bwdf_min",
                "bwdf_max",
                "bwdf_skew",
                "bwdf_kurtosis",
                "bwdf_w_mean",
                "bwdf_w_std",
            ]
            assert sorted(df_icoop.columns) == sorted(stats_df_expected_columns)

        elif bwdf_df_type == "sorted_bwdf":
            assert not df_icoop.isna().to_numpy().any()
            unique_values, counts = unique(df_icoop.values, return_counts=True)
            assert dict(zip(unique_values, counts))[0.0] == 3
            assert set(df_icoop.columns) == {f"bwdf_at_dist{c_idx}" for c_idx in range(4)}

        else:  # sorted_dists
            assert not df_icoop.isna().to_numpy().any()
            non_zero = count_nonzero(df_icoop.values)
            if sorted_dists_mode == "positive":
                assert non_zero == 4
                assert set(df_icoop.columns) == {f"dist_at_{sorted_dists_mode[:3]}_bwdf{c_idx}" for c_idx in range(2)}
            else:
                assert non_zero == 5
                assert set(df_icoop.columns) == {f"dist_at_{sorted_dists_mode[:3]}_bwdf{c_idx}" for c_idx in range(2)}


class TestExceptions:
    def test_batch_summary_featurizer_exception(self):
        with pytest.raises(ValueError) as err0:  # noqa: PT012, PT011
            self.summary_featurize_with_json_ex = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs_exceptions/1/",
                bonds="all",
                feature_type="nonbonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json_ex.get_df()

        assert str(err0.value) == (
            "Parameter feature_type set to nonbonding but must be in ['bonding', 'antibonding', 'overall']."
        )

        with pytest.raises(Exception) as err1:  # noqa: PT012, PT011
            self.summary_featurize_with_json_ex = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs_exceptions/1/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json_ex.get_df()

        assert str(err1.value) == "Files ['COBICAR.lobster', 'ICOBILIST.lobster'] not found in mp-2176."

        with pytest.raises(Exception) as err2:  # noqa: PT012, PT011
            self.summary_featurize_with_json_ex2 = BatchSummaryFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs_exceptions/2/",
                bonds="all",
                feature_type="antibonding",
                include_cobi_data=True,
                include_coop_data=True,
                e_range=[-15, 0],
            )

            _ = self.summary_featurize_with_json_ex2.get_df()

        assert str(err2.value) == "Files ['COOPCAR.lobster', 'ICOOPLIST.lobster'] not found in mp-1000."

        # COXX exception
        with pytest.raises(Exception) as err3:  # noqa: PT012, PT011
            self.raise_coxx_exception = BatchSummaryFeaturizer(path_to_lobster_calcs=TestDir / "test_data/JSONS/")

            _ = self.raise_coxx_exception._featurizecoxx(
                path_to_lobster_calc=self.raise_coxx_exception.path_to_lobster_calcs
            )

        assert str(err3.value) == "Files ['CONTCAR', 'COHPCAR.lobster', 'ICOHPLIST.lobster'] not found in JSONS."

        # Charges exception
        with pytest.raises(Exception) as err4:  # noqa: PT012, PT011
            self.raise_ch_exception = BatchSummaryFeaturizer(path_to_lobster_calcs=TestDir / "test_data/JSONS/")

            _ = self.raise_ch_exception._featurizecharges(
                path_to_lobster_calc=self.raise_ch_exception.path_to_lobster_calcs
            )

        assert str(err4.value) == "Files ['CONTCAR', 'CHARGE.lobster'] not found in JSONS."

        # Fingerprint similarity exception
        with pytest.raises(Exception) as err8:  # noqa: PT012, PT011
            fp_cohp_bonding = BatchCoxxFingerprint(
                path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
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
            str(err8.value) == "Cannot compute similarity index. Please set either normalize=True or "
            "tanimoto=True or both to False."
        )

        with pytest.raises(ValueError, match=re.escape("param mode must be in ('positive', 'negative')")):  # noqa: PT012
            batch_icoop = BatchIcoxxlistFeaturizer(
                path_to_lobster_calcs=TestDir / "test_data/Featurizer_test_data/Lobster_calcs",
                n_jobs=3,
                bin_width=0.5,
                read_icoops=True,
                normalization="area",
                bwdf_df_type="sorted_dists",
                sorted_dists_mode="invalid",
            )
            _ = batch_icoop.get_df()
