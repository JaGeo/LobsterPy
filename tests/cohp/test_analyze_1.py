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
        assert self.analyse_NaCl.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
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
        assert self.analyse_NaCl_valences.condensed_bonding_analysis["sites"][0]["bonds"]["Cl"][
            "has_antibdg_states_below_Efermi"
        ]
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
