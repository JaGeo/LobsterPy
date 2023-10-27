import pytest
from pathlib import Path

from lobsterpy.structuregraph.graph import LobsterGraph

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestGraph:
    def setup_method(self):
        self.graph_NaCl_all = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=None,
        )

        self.graph_NaCl_cation_anion = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="cation-anion",
            start=None,
        )

        self.graph_NaCl_without_add_data = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=None,
        )

        self.graph_NaCl_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=-4,
        )

        self.graph_CdF_all = LobsterGraph(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "test_data/CdF_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/CdF_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=None,
        )

        self.graph_CdF_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "test_data/CdF_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/CdF_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=-2,
        )

    def test_graph_na_cl_all(self):
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["coords"].tolist() == [2.845847, 2.845847, 2.845847]
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08484
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.02826
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOBI"] == 0.08482
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOOP"] == 0.02824
        assert self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_na_cl_cation_anion(self):
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert "env" not in self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["properties"]

        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(
            -0.5661, abs=0.001
        )
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08484
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.02826
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(
            -0.56614, abs=0.001
        )
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOBI"] == 0.08482
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOOP"] == 0.02824
        assert self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_na_cl_without_add_data(self):
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(
            -0.5661, abs=0.001
        )
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        assert "ICOOP" not in self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(
            -0.5661, abs=0.001
        )
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        assert "ICOOP" not in self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        assert self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_na_cl_close_fermi(self):
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        assert "ICOOP" not in self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        assert "ICOOP" not in self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        assert self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_cd_f_all(self):
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[0]["specie"] == "Cd"
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 1.57
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[0]["properties"]["env"] == "C:8"

        assert self.graph_CdF_all.sg.graph.nodes(data=True)[1]["specie"] == "F"
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[1]["coords"].tolist() == [1.37314, 1.37314, 1.37314]
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.79
        assert self.graph_CdF_all.sg.graph.nodes(data=True)[1]["properties"]["env"] == "T:4"

        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == -0.62168
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 0.72263
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0.27737
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08932
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.0148
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "29"

        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP"] == -0.62168
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP_bonding_perc"] == 0.72263
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP_antibonding_perc"] == 0.27737
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOBI"] == 0.08932
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOOP"] == 0.0148
        assert self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["bond_label"] == "63"

    def test_graph_cd_f_close_fermi(self):
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["specie"] == "Cd"
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 1.57
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["env"] == "C:8"

        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["specie"] == "F"
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist() == [1.37314, 1.37314, 1.37314]
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.79
        assert self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["env"] == "T:4"

        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == -0.62168
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 0
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 1
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08932
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.0148
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "29"

        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP"] == -0.62168
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP_bonding_perc"] == 0
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP_antibonding_perc"] == 1
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOBI"] == 0.08932
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOOP"] == 0.0148
        assert self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["bond_label"] == "63"
