from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest
from pymatgen.core import Structure
from pymatgen.io.lobster import Bandoverlaps, Charge, Doscar, Lobsterin, Lobsterout
from pymatgen.io.vasp import Vasprun

from lobsterpy.cohp.analyze import Analysis

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestAnalyse:
    def test_all_attributes_nacl_mulliken(self, analyse_nacl):
        assert analyse_nacl.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert analyse_nacl.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(analyse_nacl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]) == pytest.approx(
            -0.57
        )
        assert float(analyse_nacl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]) == pytest.approx(
            -3.39
        )
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["has_antibdg_states_below_Efermi"]
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert analyse_nacl.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert analyse_nacl.condensed_bonding_analysis["type_charges"] == "Mulliken"

        # capture and test printed text when orbital resolved analysis is not switched on
        captured_output = io.StringIO()
        sys.stdout = captured_output

        analyse_nacl.get_site_orbital_resolved_labels()

        console_text_printed = captured_output.getvalue().strip()

        ref_text = "Please set orbital_resolved to True when instantiating Analysis object, to get this data"
        assert console_text_printed == ref_text

    def test_all_attributes_nacl_valences(self, analyse_nacl_valences):
        assert analyse_nacl_valences.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_valences.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert analyse_nacl_valences.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]
        ) == pytest.approx(-0.57)
        assert float(
            analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1)
        assert analyse_nacl_valences.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert analyse_nacl_valences.condensed_bonding_analysis["type_charges"] == "Valences"

    def test_all_attributes_batio3(self, analyse_bati03):
        assert analyse_bati03.condensed_bonding_analysis["formula"] == "BaTiO3"
        assert analyse_bati03.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.94882)
        assert analyse_bati03.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(analyse_bati03.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["ICOHP_sum"]) == pytest.approx(
            -21.24
        )
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["has_antibdg_states_below_Efermi"]
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["bonds"]["O"]["number_of_bonds"] == 6
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["ion"] == "Ti"
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.96)
        assert analyse_bati03.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "87",
            "88",
            "98",
            "101",
            "109",
            "114",
        ]
        assert analyse_bati03.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_nacl_all(self, analyse_nacl_all):
        assert analyse_nacl_all.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_all.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_all.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert analyse_nacl_all.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            analyse_nacl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_mean"]
        ) == pytest.approx(-0.57)
        assert float(
            analyse_nacl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["has_antibdg_states_below_Efermi"]
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert analyse_nacl_all.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert analyse_nacl_all.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_nacl_comp_range_cobi_orbital(self, analyse_nacl_comp_range_cobi_orb):
        assert analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis[
            "number_of_considered_ions"
        ] == pytest.approx(1)
        assert (
            analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["ICOBI_mean"]
            == 0.0314
        )
        assert float(
            analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3s-3s"
            ]["bonding"]["integral"]
        ) == pytest.approx(0.21)
        assert float(
            analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3p-3s"
            ]["ICOBI_sum"]
        ) == pytest.approx(0.3206)
        assert float(
            analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                "3p-3s"
            ]["orb_contribution_perc_bonding"]
        ) == pytest.approx(0.6)
        assert analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
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
                analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                    "3p-3s"
                ]["bonding"]["perc"]
            )
            == 1.0
        )
        assert analyse_nacl_comp_range_cobi_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
            "3p-3s"
        ]["relevant_sub_orbitals"] == ["3py-3s", "3pz-3s", "3px-3s"]

    def test_all_attributes_nacl_comp_range_orbital(self, analyse_nacl_comp_range_orb):
        assert analyse_nacl_comp_range_orb.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert (
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3s-3s"][
                "ICOHP_mean"
            ]
            == -0.3249
        )
        assert float(
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3s-3s"][
                "bonding"
            ]["integral"]
        ) == pytest.approx(2.18)
        assert float(
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3s-3s"][
                "antibonding"
            ]["integral"]
        ) == pytest.approx(0.23)
        assert float(
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3p-3s"][
                "ICOHP_sum"
            ]
        ) == pytest.approx(-1.4485)
        assert float(
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3p-3s"][
                "orb_contribution_perc_bonding"
            ]
        ) == pytest.approx(0.39)
        assert float(
            analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"]["3p-2p"][
                "orb_contribution_perc_antibonding"
            ]
        ) == pytest.approx(0.32)
        assert analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
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
                analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                    "3p-3s"
                ]["bonding"]["perc"]
            )
            == 0.9932
        )
        assert (
            float(
                analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
                    "3p-2p"
                ]["antibonding"]["perc"]
            )
            == 0.5
        )
        assert analyse_nacl_comp_range_orb.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["orbital_data"][
            "3p-2p"
        ]["relevant_sub_orbitals"] == [
            "3py-2py",
            "3pz-2py",
            "3px-2py",
            "3py-2pz",
            "3pz-2pz",
            "3px-2pz",
            "3py-2px",
            "3pz-2px",
            "3px-2px",
        ]

    def test_all_attributes_analyse_nacl_comp_range_cobi(self, analyse_nacl_comp_range_cobi):
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            5.69169
        )
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOBI_mean"]
        ) == pytest.approx(0.08)
        assert float(
            analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOBI_sum"]
        ) == pytest.approx(0.51)
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
        assert (
            analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        )
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert analyse_nacl_comp_range_cobi.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_final_dicts(
        self,
        analyse_nasi_madelung_all,
        analyse_nacl,
        analyse_bati03,
        analyse_nacl_distorted,
        analyse_nacl_spin,
        analyse_nacl_all,
    ):
        assert analyse_nasi_madelung_all.final_dict_bonds["Na-Si"]["ICOHP_mean"] == pytest.approx(-0.44833333333)
        assert analyse_nacl.final_dict_ions["Na"] == {"O:6": 1}
        assert analyse_nacl.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert analyse_nacl.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.565000000000000)

        assert analyse_bati03.final_dict_ions["Ti"] == {"O:6": 1}
        assert analyse_bati03.final_dict_bonds["O-Ti"]["has_antbdg"]
        assert analyse_bati03.final_dict_bonds["O-Ti"]["ICOHP_mean"] == pytest.approx(-3.5399999999999996)

        assert analyse_nacl_distorted.final_dict_ions == {"Na": {"O:6": 4}}
        assert analyse_nacl_distorted.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert analyse_nacl_distorted.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5658333333333334)
        assert analyse_nacl_spin.final_dict_bonds["Cl-Na"]["has_antbdg"]
        assert analyse_nacl_spin.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5650000000000001)
        assert analyse_nacl_spin.final_dict_ions == {"Na": {"O:6": 1}}
        assert analyse_nacl_all.final_dict_bonds["Cl-Na"]["ICOHP_mean"] == pytest.approx(-0.5650000000000001)
        assert analyse_nacl_all.final_dict_ions == {
            "Na": {"O:6": 1},
            "Cl": {"O:6": 1},
        }

    def test_all_attributes_nasbf6(self, analyse_nasbf6):
        assert analyse_nasbf6.condensed_bonding_analysis["formula"] == "NaSbF6"
        assert analyse_nasbf6.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.98893)
        assert analyse_nasbf6.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]) == pytest.approx(
            -3.65
        )
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 6
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["integral"] == 3.77
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["perc"] == 0.95929
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 0.04071
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["integral"] == 0.16
        assert abs(
            float(analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["integral"]
                    - analyse_nasbf6.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["integral"],
                    2,
                )
            ),
            abs=0.10,
        )
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.91)
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "25",
            "31",
            "34",
            "43",
            "47",
        ]
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"]) == pytest.approx(
            -32.71
        )
        assert not (
            analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        )
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["number_of_bonds"] == 6
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"]["perc"] == pytest.approx(
            1.0
        )
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        assert abs(
            float(analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"]["integral"]
                    - analyse_nasbf6.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"]["integral"],
                    2,
                )
            ),
            abs=0.10,
        )
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["ion"] == "Sb"
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(2.91)
        assert analyse_nasbf6.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "63",
            "69",
            "73",
            "80",
            "81",
            "87",
        ]
        assert analyse_nasbf6.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_nasbf6_anbd(self, analyse_nasbf6_anbd):
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["formula"] == "NaSbF6"
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.98893)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-3.65)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 6
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(0.0)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.91)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "25",
            "31",
            "34",
            "43",
            "47",
        ]
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["env"] == "O:6"
        assert float(
            analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["ICOHP_sum"]
        ) == pytest.approx(-32.71)
        assert not (
            analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        )
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["number_of_bonds"] == 6
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(1.0)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["bonds"]["F"]["antibonding"]["perc"] == 0.0
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["ion"] == "Sb"
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(2.91)
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "63",
            "69",
            "73",
            "80",
            "81",
            "87",
        ]
        assert analyse_nasbf6_anbd.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_cdf(self, analyse_cdf):
        assert analyse_cdf.condensed_bonding_analysis["formula"] == "CdF2"
        assert analyse_cdf.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.99501)
        assert analyse_cdf.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["env"] == "C:8"
        assert float(analyse_cdf.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"]) == pytest.approx(
            -4.96
        )
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["has_antibdg_states_below_Efermi"]
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 8
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["perc"] == pytest.approx(0.0)
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"] == 1.0
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["ion"] == "Cd"
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1.57)
        assert analyse_cdf.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "25",
            "32",
            "35",
            "36",
            "57",
            "58",
            "61",
            "68",
        ]
        assert analyse_cdf.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_cdf_comp_range_coop(self, analyse_cdf_comp_range_coop):
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["formula"] == "CdF2"
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            5.98538
        )
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["env"] == "C:8"
        assert float(
            analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOOP_sum"]
        ) == pytest.approx(0.12)
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["number_of_bonds"] == 8
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"][
            "perc"
        ] == pytest.approx(0.59016)
        assert (
            analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"]["perc"]
            == 0.40984
        )
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["ion"] == "Cd"
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(1.57)
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "29",
            "30",
            "33",
            "40",
            "53",
            "60",
            "63",
            "64",
        ]
        assert analyse_cdf_comp_range_coop.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_k3sb(self, analyse_k3sb):
        assert analyse_k3sb.condensed_bonding_analysis["formula"] == "K3Sb"
        assert analyse_k3sb.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(4.28164)
        assert analyse_k3sb.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["env"] == "6"
        assert float(analyse_k3sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOHP_sum"]) == pytest.approx(
            -0.84
        )
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert analyse_k3sb.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
        ]

        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["env"] == "4"
        assert float(analyse_k3sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOHP_sum"]) == pytest.approx(
            -1.45
        )
        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["number_of_bonds"] == 4
        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert analyse_k3sb.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == ["21", "22", "23", "24"]

        assert analyse_k3sb.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_k3sb_all(self, analyse_k3sb_all):
        assert analyse_k3sb_all.condensed_bonding_analysis["formula"] == "K3Sb"
        assert analyse_k3sb_all.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(4.28164)
        assert analyse_k3sb_all.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(3)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["env"] == "14"
        assert float(
            analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-2.97)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["number_of_bonds"] == 8
        assert float(
            analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-0.84)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
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

        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["env"] == "14"
        assert float(
            analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOHP_sum"]
        ) == pytest.approx(-1.45)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 10
        assert float(
            analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-2.15)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 10
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
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

        assert analyse_k3sb_all.condensed_bonding_analysis["type_charges"] == "Mulliken"

        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["env"] == "14"
        assert float(
            analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOHP_sum"]
        ) == pytest.approx(-3.74)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["has_antibdg_states_below_Efermi"]
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 14
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert analyse_k3sb_all.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
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

    def test_all_attributes_k3sb_all_cobi(self, analyse_k3sb_all_cobi):
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["formula"] == "K3Sb"
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(4.28164)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(3)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["env"] == "14"
        assert float(
            analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(0.26)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["K"]["number_of_bonds"] == 8
        assert float(
            analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["ICOBI_sum"]
        ) == pytest.approx(0.41)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["bonds"]["Sb"]["number_of_bonds"] == 6
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["ion"] == "K"
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.68)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
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

        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["env"] == "8"
        assert float(
            analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOBI_sum"]
        ) == pytest.approx(0.54)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 4
        assert float(
            analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(0.13)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["bonds"]["K"]["number_of_bonds"] == 4
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "1",
            "2",
            "3",
            "4",
            "21",
            "22",
            "23",
            "24",
        ]

        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["type_charges"] == "Mulliken"

        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["env"] == "14"
        assert float(
            analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOBI_sum"]
        ) == pytest.approx(1.48)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 14
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert analyse_k3sb_all_cobi.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
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

    def test_all_attributes_k3sb_all_coop_orb(self, analyse_k3sb_all_coop_orb):
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["formula"] == "K3Sb"
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(
            4.28164
        )
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(2)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["env"] == "T:4"
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["ICOOP_sum"]
        ) == pytest.approx(0.29)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["ion"] == "K"
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["charge"] == pytest.approx(0.52)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["relevant_bonds"] == [
            "21",
            "22",
            "23",
            "24",
        ]

        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["type_charges"] == "Mulliken"
        assert (
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"]["5s-4s"][
                "ICOOP_mean"
            ]
            == 0.0251
        )
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"]["5s-4s"][
                "bonding"
            ]["integral"]
        ) == pytest.approx(0.12)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"]["5s-4s"][
                "antibonding"
            ]["integral"]
        ) == pytest.approx(0.02)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"]["5p-4s"][
                "ICOOP_sum"
            ]
        ) == pytest.approx(0.2383)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"]["5p-4s"][
                "orb_contribution_perc_bonding"
            ]
        ) == pytest.approx(0.65)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
            "relevant_bonds"
        ] == ["21", "22", "23", "24"]
        assert (
            float(
                analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                    "5p-4s"
                ]["bonding"]["perc"]
            )
            == 0.96
        )
        assert (
            float(
                analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][1]["bonds"]["Sb"]["orbital_data"][
                    "5s-4s"
                ]["antibonding"]["perc"]
            )
            == 0.14286
        )

        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["env"] == "C:8"
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["ICOOP_sum"]
        ) == pytest.approx(0.59)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"][
            "has_antibdg_states_below_Efermi"
        ]
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["number_of_bonds"] == 8
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["ion"] == "Sb"
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["charge"] == pytest.approx(-1.73)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"]["5p-4s"][
                "orb_contribution_perc_bonding"
            ]
        ) == pytest.approx(0.66)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"]["5p-4s"][
                "bonding"
            ]["integral"]
        ) == pytest.approx(0.49)
        assert float(
            analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["bonds"]["K"]["orbital_data"]["5s-4s"][
                "bonding"
            ]["perc"]
        ) == pytest.approx(0.88889)
        assert analyse_k3sb_all_coop_orb.condensed_bonding_analysis["sites"][3]["relevant_bonds"] == [
            "21",
            "22",
            "23",
            "24",
            "25",
            "26",
            "27",
            "28",
        ]

    def test_all_attributes_nacl_nan(self, analyse_nacl_nan):
        assert analyse_nacl_nan.condensed_bonding_analysis["formula"] == "NaCl"
        assert analyse_nacl_nan.condensed_bonding_analysis["max_considered_bond_length"] == pytest.approx(5.69169)
        assert analyse_nacl_nan.condensed_bonding_analysis["number_of_considered_ions"] == pytest.approx(1)
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["env"] == "O:6"
        assert float(
            analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"]
        ) == pytest.approx(-3.39)
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["has_antibdg_states_below_Efermi"]
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["number_of_bonds"] == 6
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["bonding"][
            "perc"
        ] == pytest.approx(0.0)
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["antibonding"]["perc"] == 0.0
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["ion"] == "Na"
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["charge"] == pytest.approx(0.78)
        assert analyse_nacl_nan.condensed_bonding_analysis["sites"][0]["relevant_bonds"] == [
            "21",
            "23",
            "24",
            "27",
            "28",
            "30",
        ]
        assert analyse_nacl_nan.condensed_bonding_analysis["type_charges"] == "Mulliken"

    def test_all_attributes_nacl_pymatgen_objs(
        self,
        analyse_nacl_comp_range_orb,
        analyse_nacl_comp_range_orb_with_objs,
        analyse_k3sb_all,
        analyse_k3sb_all_objs,
    ):
        assert (
            analyse_nacl_comp_range_orb.condensed_bonding_analysis
            == analyse_nacl_comp_range_orb_with_objs.condensed_bonding_analysis
        )
        assert analyse_k3sb_all.condensed_bonding_analysis == analyse_k3sb_all_objs.condensed_bonding_analysis

    def test_icohp_sum_nacl(self, analyse_nacl_comp_range):
        assert abs(
            float(analyse_nacl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    analyse_nacl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["bonding"]["integral"]
                    - analyse_nacl_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )

    def test_icohp_sum_cdf(self, analyse_cdf_comp_range):
        assert abs(
            float(analyse_cdf_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["ICOHP_sum"])
        ) == pytest.approx(
            (
                round(
                    analyse_cdf_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["bonding"]["integral"]
                    - analyse_cdf_comp_range.condensed_bonding_analysis["sites"][0]["bonds"]["F"]["antibonding"][
                        "integral"
                    ],
                    2,
                )
            ),
            abs=0.10,
        )

    def test_exception(self):
        with pytest.raises(ValueError):  # noqa: PT011
            self.analyse_batao2n1 = Analysis(
                path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
                path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
                path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
                which_bonds="cation_cation",
                cutoff_icohp=0.1,
            )
        with pytest.raises(ValueError) as err:  # noqa: PT011
            self.analyse_c = Analysis(
                path_to_poscar=TestDir / "test_data/C/POSCAR.gz",
                path_to_cohpcar=TestDir / "test_data/C/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "test_data/C/ICOHPLIST.lobster.gz",
                path_to_charge=TestDir / "test_data/C/CHARGE.lobster.gz",
                which_bonds="cation-anion",
                cutoff_icohp=0.1,
            )
        assert (
            str(err.value) == "Consider switching to an analysis of all bonds and not only cation-anion bonds. "
            "It looks like no cations are detected."
        )


class TestAnalyseCalcQuality:
    def test_calc_quality_summary_with_objs(self):
        charge_obj = Charge(filename=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz")
        bandoverlaps_obj = Bandoverlaps(filename=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz")
        structure_obj = Structure.from_file(filename=TestDir / "test_data" / "K3Sb" / "POSCAR.gz")
        vasprun_obj = Vasprun(
            filename=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz", parse_eigen=False, parse_potcar_file=False
        )
        doscar = Doscar(
            doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            structure_file=None,
            structure=structure_obj,
        )
        lobsterin_obj = Lobsterin.from_file(TestDir / "test_data" / "K3Sb" / "lobsterin.gz")
        lobsterout_obj = Lobsterout(filename=TestDir / "test_data" / "K3Sb" / "lobsterout.gz")

        calc_des_with_objs = Analysis.get_lobster_calc_quality_summary(
            structure_obj=structure_obj,
            lobster_completedos_obj=doscar.completedos,
            charge_obj=charge_obj,
            bandoverlaps_obj=bandoverlaps_obj,
            vasprun_obj=vasprun_obj,
            lobsterin_obj=lobsterin_obj,
            lobsterout_obj=lobsterout_obj,
            dos_comparison=True,
            bva_comp=True,
            n_bins=256,
        )

        calc_des_with_paths = Analysis.get_lobster_calc_quality_summary(
            path_to_poscar=TestDir / "test_data" / "K3Sb" / "POSCAR.gz",
            potcar_symbols=["K_sv", "Sb"],
            path_to_charge=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz",
            path_to_doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            path_to_vasprun=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz",
            path_to_bandoverlaps=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz",
            path_to_lobsterout=TestDir / "test_data" / "K3Sb" / "lobsterout.gz",
            path_to_lobsterin=TestDir / "test_data" / "K3Sb" / "lobsterin.gz",
            dos_comparison=True,
            bva_comp=True,
            n_bins=256,
        )

        assert calc_des_with_objs == calc_des_with_paths

    def test_calc_quality_summary_exceptions(self):
        charge_obj = Charge(filename=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz")
        bandoverlaps_obj = Bandoverlaps(filename=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz")
        structure_obj = Structure.from_file(filename=TestDir / "test_data" / "K3Sb" / "POSCAR.gz")
        vasprun_obj = Vasprun(
            filename=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz", parse_eigen=False, parse_potcar_file=False
        )
        doscar = Doscar(
            doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            structure_file=None,
            structure=structure_obj,
        )
        lobsterin_obj = Lobsterin.from_file(TestDir / "test_data" / "K3Sb" / "lobsterin.gz")
        lobsterout_obj = Lobsterout(filename=TestDir / "test_data" / "K3Sb" / "lobsterout.gz")

        with pytest.raises(ValueError) as err_poscar:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=None,
                structure_obj=None,
                lobster_completedos_obj=doscar.completedos,
                charge_obj=charge_obj,
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=vasprun_obj,
                lobsterin_obj=lobsterin_obj,
                lobsterout_obj=lobsterout_obj,
            )

        assert str(err_poscar.value) == "Please provide path_to_poscar or structure_obj"

        with pytest.raises(ValueError) as err_potcar:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=None,
                path_to_potcar=None,
                potcar_symbols=None,
                structure_obj=structure_obj,
                lobster_completedos_obj=doscar.completedos,
                charge_obj=charge_obj,
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=None,
                lobsterin_obj=lobsterin_obj,
                lobsterout_obj=lobsterout_obj,
            )

        assert (
            str(err_potcar.value) == "Please provide either path_to_potcar or list of "
            "potcar_symbols or path to vasprun.xml or vasprun object. "
            "Crucial to identify basis used for projections"
        )

        with pytest.raises(ValueError) as err_lobsterout:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                structure_obj=structure_obj,
                lobster_completedos_obj=doscar.completedos,
                charge_obj=charge_obj,
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=vasprun_obj,
                lobsterin_obj=lobsterin_obj,
                lobsterout_obj=None,
            )

        assert str(err_lobsterout.value) == "Please provide path_to_lobsterout or lobsterout_obj"

        with pytest.raises(ValueError) as err_lobsterin:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                structure_obj=structure_obj,
                lobster_completedos_obj=doscar.completedos,
                charge_obj=charge_obj,
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=vasprun_obj,
                lobsterin_obj=None,
                lobsterout_obj=lobsterout_obj,
            )

        assert str(err_lobsterin.value) == "Please provide path_to_lobsterin or lobsterin_obj"

        with pytest.raises(Exception) as err_charge:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                structure_obj=structure_obj,
                lobster_completedos_obj=doscar.completedos,
                charge_obj=None,
                path_to_charge=None,
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=vasprun_obj,
                lobsterin_obj=lobsterin_obj,
                lobsterout_obj=lobsterout_obj,
                bva_comp=True,
            )

        assert str(err_charge.value) == "BVA comparison is requested, thus please provide path_to_charge or charge_obj"

        with pytest.raises(ValueError) as err_doscar:  # noqa: PT011
            self.calc_des = Analysis.get_lobster_calc_quality_summary(
                structure_obj=structure_obj,
                lobster_completedos_obj=None,
                charge_obj=charge_obj,
                potcar_symbols=["K_sv", "Sb"],
                bandoverlaps_obj=bandoverlaps_obj,
                vasprun_obj=None,
                lobsterin_obj=lobsterin_obj,
                lobsterout_obj=lobsterout_obj,
                bva_comp=True,
                dos_comparison=True,
            )

        assert (
            str(err_doscar.value)
            == "Dos comparison is requested, so please provide either path_to_doscar or lobster_completedos_obj"
        )
