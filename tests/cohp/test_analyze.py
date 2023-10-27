from __future__ import annotations

import pytest
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestAnalyse:
    def setup_method(self):
        self.analyse_NaCl = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_comp_range = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_comp_range_orb = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            orbital_cutoff=0.10,
            orbital_resolved=True,
        )
        self.analyse_NaCl_comp_range_cobi = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            noise_cutoff=0.001,
            are_cobis=True,
        )
        self.analyse_NaCl_comp_range_cobi_orb = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            noise_cutoff=0.001,
            are_cobis=True,
            orbital_resolved=True,
        )

        self.analyse_NaCl_nan = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
        )

        self.analyse_NaCl_valences = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=None,
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_madelung = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            path_to_madelung=TestDir / "test_data/NaCl/MadelungEnergies.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTiO3 = Analysis(
            path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTaO2N1 = Analysis(
            path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTiO3_differentcutoff = Analysis(
            path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.001,
        )

        self.analyse_NaCl_distorted = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_distorted/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl_distorted/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl_distorted/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl_distorted/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_spin = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_spin/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl_spin/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl_spin/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl_spin/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_all = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            which_bonds="all",
            cutoff_icohp=0.1,
        )

        self.analyse_NaCl_madelung_all = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            path_to_madelung=TestDir / "test_data/NaCl/MadelungEnergies.lobster",
            which_bonds="all",
            cutoff_icohp=0.1,
        )
        self.analyse_NaSi_madelung_all = Analysis(
            path_to_poscar=TestDir / "test_data/NaSi/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaSi/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaSi/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaSi/CHARGE.lobster",
            path_to_madelung=TestDir / "test_data/NaSi/MadelungEnergies.lobster",
            which_bonds="all",
            cutoff_icohp=0.1,
        )

        self.analyse_BaTaO2N1_cutoff = Analysis(
            path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.001,
        )

        self.analyse_NaSbF6 = Analysis(
            path_to_poscar=TestDir / "test_data/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_NaSbF6_anbd = Analysis(
            path_to_poscar=TestDir / "test_data/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            start=-5.5,
        )

        self.analyse_CdF = Analysis(
            path_to_poscar=TestDir / "test_data/CdF/POSCAR",
            path_to_cohpcar=TestDir / "test_data/CdF/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/CdF/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
        )

        self.analyse_CdF_comp_range = Analysis(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_CdF_comp_range_coop = Analysis(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COOPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            noise_cutoff=0.001,
            are_coops=True,
        )

        self.analyse_K3Sb = Analysis(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.analyse_K3Sb_all = Analysis(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
        )

        self.analyse_K3Sb_all_cobi = Analysis(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/K3Sb/COBICAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/K3Sb/ICOBILIST.lobster.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
            noise_cutoff=0.001,
            are_cobis=True,
        )

        self.analyse_K3Sb_all_coop_orb = Analysis(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/K3Sb/COOPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/K3Sb/ICOOPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
            noise_cutoff=0.001,
            orbital_resolved=True,
            are_coops=True,
        )

        # different environment than O:6

    def test_exception(self):
        with pytest.raises(ValueError):
            self.analyse_BaTaO2N1 = Analysis(
                path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
                path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
                path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
                which_bonds="cation_cation",
                cutoff_icohp=0.1,
            )
        with pytest.raises(ValueError) as err:
            self.analyse_C = Analysis(
                path_to_poscar=TestDir / "test_data/C/POSCAR",
                path_to_cohpcar=TestDir / "test_data/C/COHPCAR.lobster",
                path_to_icohplist=TestDir / "test_data/C/ICOHPLIST.lobster",
                path_to_charge=TestDir / "test_data/C/CHARGE.lobster",
                which_bonds="cation-anion",
                cutoff_icohp=0.1,
            )
        assert (
            str(err.value) == "Consider switching to an analysis of all bonds and not only cation-anion bonds. "
            "It looks like no cations are detected."
        )

    def test_all_attributes_na_cl_mulliken(self):
        assert self.analyse_NaCl.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert self.analyse_NaCl.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]
        ) == pytest.approx(-0.57)
        assert float(
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert (
            self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["has_antibdg_states_below_Efermi"]
        )
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert self.analyse_NaCl.condensed_bonding_analysis["type_charges"] == "Mulliken"
        assert self.analyse_NaCl_madelung.condensed_bonding_analysis["madelung_energy"] == pytest.approx(-5.40)

    def test_all_attributes_na_cl_valences(self):
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            5.69169
        )
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]
        ) == pytest.approx(-0.57)
        assert float(
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert (
            self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "has_antibdg_states_below_Efermi"]
        )
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1)
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["type_charges"] == "Valences"

    def test_all_attributes_ba_ti_o3(self):
        assert self.analyse_BaTiO3.condensed_bonding_analysis["formula"] == "BaTiO3"
        assert self.analyse_BaTiO3.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.94882)
        assert self.analyse_BaTiO3.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["ICOHP_sum"]
        ) == pytest.approx(-21.24)
        assert (
            self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["has_antibdg_states_below_Efermi"]
        )
        assert self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["number_of_bonds"] == 6
        assert self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["ion"] == "Ti"
        assert self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.96)
        assert self.analyse_BaTiO3.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "87",
            "88",
            "98",
            "101",
            "109",
            "114",
        ]
        assert self.analyse_BaTiO3.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_na_cl_all(self):
        assert self.analyse_NaCl_all.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_all.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_all.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert self.analyse_NaCl_all.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]
        ) == pytest.approx(-0.57)
        assert float(
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert (
            self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                "has_antibdg_states_below_Efermi"
            ]
        )
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert self.analyse_NaCl_all.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert self.analyse_NaCl_all.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_na_cl_comp_range_cobi_orbital(self):
        assert self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis[
            "number_of_considered_ions"
        ] == pytest.approx(1)
        assert (
            self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["ICOBI_mean"]
            == 0.0314
        )
        assert float(
            self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["bonding"]["integral"]
        ) == pytest.approx(0.21)
        assert float(
            self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3px-3s"
            ]["ICOBI_sum"]
        ) == pytest.approx(0.1069)
        assert float(
            self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3py-3s"
            ]["orb_contribution_perc_bonding"]
        ) == pytest.approx(0.2)
        assert self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "orbital_data"
        ]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert (
            float(
                self.analyse_NaCl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
                    "orbital_data"
                ]["3pz-3s"]["bonding"]["perc"]
            )
            == 1.0
        )

    def test_all_attributes_na_cl_comp_range_orbital(self):
        assert self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis[
            "number_of_considered_ions"
        ] == pytest.approx(1)
        assert (
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["ICOHP_mean"]
            == -0.3249
        )
        assert float(
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["bonding"]["integral"]
        ) == pytest.approx(2.18)
        assert float(
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["antibonding"]["integral"]
        ) == pytest.approx(0.23)
        assert float(
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3px-3s"
            ]["ICOHP_sum"]
        ) == pytest.approx(-0.4828)
        assert float(
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3py-3s"
            ]["orb_contribution_perc_bonding"]
        ) == pytest.approx(0.13)
        assert float(
            self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3px-2px"
            ]["orb_contribution_perc_antibonding"]
        ) == pytest.approx(0.11)
        assert self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
            "relevant_bonds"
        ] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert (
            float(
                self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                    "3pz-3s"
                ]["bonding"]["perc"]
            )
            == 1.0
        )
        assert (
            float(
                self.analyse_NaCl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                    "3pz-2pz"
                ]["antibonding"]["perc"]
            )
            == 0.5
        )

    def test_all_attributes_analyse_na_cl_comp_range_cobi(self):
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis[
            "max_considered_bond_length"
        ] == pytest.approx(5.69169)
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis[
            "number_of_considered_ions"
        ] == pytest.approx(1)
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOBI_mean"]
        ) == pytest.approx(0.08)
        assert float(
            self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOBI_sum"]
        ) == pytest.approx(0.51)
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
        assert (
            self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"]
            == 6
        )
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert self.analyse_NaCl_comp_range_cobi.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_final_dicts(self):
        assert self.analyse_NaSi_madelung_all.final_dict_bonds["Na-Si"]["ICOHP_mean"] == pytest.approx(-0.44833333333)
        assert self.analyse_NaCl.final_dict_ions["Na"] == {"O:6": 1}
        assert self.analyse_NaCl.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert self.analyse_NaCl.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.565000000000000)

        assert self.analyse_BaTiO3.final_dict_ions["Ti"] == {"O:6": 1}
        assert self.analyse_BaTiO3.final_dict_bonds["O-Ti"]["has_antbdg"]
        assert self.analyse_BaTiO3.final_dict_bonds["O-Ti"]["ICOHP_mean"] == pytest.approx(-3.5399999999999996)

        assert self.analyse_NaCl_distorted.final_dict_ions == {"Na": {"O:6": 4}}
        assert self.analyse_NaCl_distorted.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert self.analyse_NaCl_distorted.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5658333333333334)
        assert self.analyse_NaCl_spin.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert self.analyse_NaCl_spin.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5650000000000001)
        assert self.analyse_NaCl_spin.final_dict_ions == {"Na": {"O:6": 1}}
        assert self.analyse_NaCl_all.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5650000000000001)
        assert self.analyse_NaCl_all.final_dict_ions == {
            "Na": {"O:6": 1},
            "Cl": {"O:6": 1},
        }

    def test_all_attributes_na_sb_f6(self):
        assert self.analyse_NaSbF6.condensed_bonding_analysis["formula"] == "NaSbF6"
        assert self.analyse_NaSbF6.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.98893)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-3.65)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 6
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["integral"] == 3.77
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["perc"] == 0.95929
        assert (
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 0.04071
        )
        assert (
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["integral"] == 0.16
        )
        assert abs(
            float(self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["integral"]
                    - self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.91)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "25",
            "31",
            "34",
            "43",
            "47",
        ]
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-32.71)
        assert not (
            self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        )
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["number_of_bonds"] == 6
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(1.0)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        assert abs(
            float(self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"]["integral"]
                    - self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["ion"] == "Sb"
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(2.91)
        assert self.analyse_NaSbF6.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "63",
            "69",
            "73",
            "80",
            "81",
            "87",
        ]
        assert self.analyse_NaSbF6.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_na_sb_f6_anbd(self):
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["formula"] == "NaSbF6"
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            5.98893
        )
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-3.65)
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 6
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(0.0)
        assert (
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        )
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.91)
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "25",
            "31",
            "34",
            "43",
            "47",
        ]
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-32.71)
        assert not (
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"][
                "has_antibdg_states_below_Efermi"
            ]
        )
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["number_of_bonds"] == 6
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(1.0)
        assert (
            self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        )
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["ion"] == "Sb"
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(2.91)
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "63",
            "69",
            "73",
            "80",
            "81",
            "87",
        ]
        assert self.analyse_NaSbF6_anbd.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_cd_f(self):
        assert self.analyse_CdF.condensed_bonding_analysis["formula"] == "CdF2"
        assert self.analyse_CdF.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.99501)
        assert self.analyse_CdF.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["env"] == "C:8"
        assert float(
            self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-4.96)
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 8
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(0.0)
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 1.0
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["ion"] == "Cd"
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1.57)
        assert self.analyse_CdF.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "25",
            "32",
            "35",
            "36",
            "57",
            "58",
            "61",
            "68",
        ]
        assert self.analyse_CdF.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_cd_f_comp_range_coop(self):
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["formula"] == "CdF2"
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis[
            "max_considered_bond_length"
        ] == pytest.approx(5.98538)
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis[
            "number_of_considered_ions"
        ] == pytest.approx(1)
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["env"] == "C:8"
        assert float(
            self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOOP_sum"]
        ) == pytest.approx(0.12)
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
            "has_antibdg_states_below_Efermi"
        ]
        assert (
            self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"]
            == 8
        )
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(0.59016)
        assert (
            self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"]
            == 0.40984
        )
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["ion"] == "Cd"
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1.57)
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "29",
            "30",
            "33",
            "40",
            "53",
            "60",
            "63",
            "64",
        ]
        assert self.analyse_CdF_comp_range_coop.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_k3_sb(self):
        assert self.analyse_K3Sb.condensed_bonding_analysis["formula"] == "K3Sb"
        assert self.analyse_K3Sb.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(4.28164)
        assert self.analyse_K3Sb.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["env"] == "6"
        assert float(
            self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-0.84)
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
        ]

        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["env"] == "4"
        assert float(
            self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-1.45)
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["number_of_bonds"] == 4
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert self.analyse_K3Sb.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == ["21", "22", "23", "24"]

        assert self.analyse_K3Sb.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_k3_sb_all(self):
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["formula"] == "K3Sb"
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(4.28164)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(3)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["env"] == "14"
        assert float(
            self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-2.97)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["number_of_bonds"] == 8
        assert float(
            self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-0.84)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
        ]

        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["env"] == "14"
        assert float(
            self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-1.45)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 10
        assert float(
            self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-2.15)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 10
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "1",
            "2",
            "3",
            "4",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "23",
            "24",
        ]

        assert self.analyse_K3Sb_all.condensed_bonding_analysis["type_charges"] == "Mulliken"

        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["env"] == "14"
        assert float(
            self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-3.74)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 14
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert self.analyse_K3Sb_all.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
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

    def test_all_attributes_k3_sb_all_cobi(self):
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["formula"] == "K3Sb"
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            4.28164
        )
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(3)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["env"] == "14"
        assert float(
            self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(0.26)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["number_of_bonds"] == 8
        assert float(
            self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOBI_sum"]
        ) == pytest.approx(0.41)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
        ]

        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["env"] == "8"
        assert float(
            self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOBI_sum"]
        ) == pytest.approx(0.54)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 4
        assert float(
            self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(0.13)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 4
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "1",
            "2",
            "3",
            "4",
            "21",
            "22",
            "23",
            "24",
        ]

        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["type_charges"] == "Mulliken"

        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["env"] == "14"
        assert float(
            self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(1.48)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 14
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert self.analyse_K3Sb_all_cobi.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
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

    def test_all_attributes_k3_sb_all_coop_orb(self):
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["formula"] == "K3Sb"
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            4.28164
        )
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(
            2
        )
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["env"] == "T:4"
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOOP_sum"]
        ) == pytest.approx(0.29)
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "21",
            "22",
            "23",
            "24",
        ]

        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["type_charges"] == "Mulliken"
        assert (
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                "5s-4s"
            ]["ICOOP_mean"]
            == 0.0251
        )
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                "5s-4s"
            ]["bonding"]["integral"]
        ) == pytest.approx(0.12)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                "5s-4s"
            ]["antibonding"]["integral"]
        ) == pytest.approx(0.02)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                "5px-4s"
            ]["ICOOP_sum"]
        ) == pytest.approx(0.0796)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                "5pz-4s"
            ]["orb_contribution_perc_bonding"]
        ) == pytest.approx(0.22)
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
            "relevant_bonds"
        ] == ["21", "22", "23", "24"]
        assert (
            float(
                self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                    "5py-4s"
                ]["bonding"]["perc"]
            )
            == 1.0
        )
        assert (
            float(
                self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                    "5s-4s"
                ]["antibonding"]["perc"]
            )
            == 0.14286
        )

        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["env"] == "C:8"
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOOP_sum"]
        ) == pytest.approx(0.59)
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert (
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 8
        )
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"][
                "5pz-4s"
            ]["orb_contribution_perc_bonding"]
        ) == pytest.approx(0.22)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"][
                "5px-4s"
            ]["bonding"]["integral"]
        ) == pytest.approx(0.16)
        assert float(
            self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"][
                "5s-4s"
            ]["bonding"]["perc"]
        ) == pytest.approx(0.88889)
        assert self.analyse_K3Sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
            "21",
            "22",
            "23",
            "24",
            "25",
            "26",
            "27",
            "28",
        ]

    def test_all_attributes_na_cl_nan(self):
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["formula"] == "NaCl"
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["bonding"][
            "perc"
        ] == pytest.approx(0.0)
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["antibonding"]["perc"] == 0.0
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert self.analyse_NaCl_nan.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_icohp_sum_na_cl(self):
        assert abs(
            float(self.analyse_NaCl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    self.analyse_NaCl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["bonding"][
                        "integral"
                    ]
                    - self.analyse_NaCl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )

    def test_icohp_sum_cd_f(self):
        assert abs(
            float(self.analyse_CdF_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    self.analyse_CdF_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
                        "integral"
                    ]
                    - self.analyse_CdF_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )
