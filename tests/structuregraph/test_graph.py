from pathlib import Path

import pytest

from lobsterpy.structuregraph.graph import LobsterGraph

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestGraph:
    def test_graph_nacl_all(self):
        graph_nacl_all = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
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
        assert graph_nacl_all.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert graph_nacl_all.sg.graph.nodes(data=True)[0]["coords"].tolist() == [
            0.0,
            0.0,
            0.0,
        ]
        assert graph_nacl_all.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert graph_nacl_all.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert graph_nacl_all.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert graph_nacl_all.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert graph_nacl_all.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert graph_nacl_all.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08484
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.02826
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["ICOBI"] == 0.08482
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["ICOOP"] == 0.02824
        assert graph_nacl_all.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_nacl_cation_anion(self):
        graph_nacl_cation_anion = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
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
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert graph_nacl_cation_anion.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert "env" not in graph_nacl_cation_anion.sg.graph.nodes(data=True)[1]["properties"]

        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08484
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.02826
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.56614, abs=0.001)
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOBI"] == 0.08482
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOOP"] == 0.02824
        assert graph_nacl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_nacl_without_add_data(self):
        graph_nacl_without_add_data = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=None,
        )
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert graph_nacl_without_add_data.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        assert "ICOOP" not in graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        assert "ICOOP" not in graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        assert graph_nacl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_nacl_close_fermi(self):
        graph_nacl_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=-4,
        )
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[0]["specie"] == "Na"
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 0.78
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["env"] == "O:6"

        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[1]["specie"] == "Cl"
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            2.845847,
            2.845847,
            2.845847,
        ]
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.78
        assert graph_nacl_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["env"] == "O:6"

        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        assert "ICOOP" not in graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"

        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP"] == pytest.approx(-0.5661, abs=0.001)
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"] == 1
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP_antibonding_perc"] == 0
        assert "ICOBI" not in graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        assert "ICOOP" not in graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        assert graph_nacl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["bond_label"] == "28"

    def test_graph_cdf_all(self):
        graph_cdf_all = LobsterGraph(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
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
        assert graph_cdf_all.sg.graph.nodes(data=True)[0]["specie"] == "Cd"
        assert graph_cdf_all.sg.graph.nodes(data=True)[0]["coords"].tolist() == [
            0.0,
            0.0,
            0.0,
        ]
        assert graph_cdf_all.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 1.57
        assert graph_cdf_all.sg.graph.nodes(data=True)[0]["properties"]["env"] == "C:8"

        assert graph_cdf_all.sg.graph.nodes(data=True)[1]["specie"] == "F"
        assert graph_cdf_all.sg.graph.nodes(data=True)[1]["coords"].tolist() == [
            1.37314,
            1.37314,
            1.37314,
        ]
        assert graph_cdf_all.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.79
        assert graph_cdf_all.sg.graph.nodes(data=True)[1]["properties"]["env"] == "T:4"

        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == -0.62168
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 0.73333
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 0.26667
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08932
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.0148
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "29"

        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP"] == -0.62168
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP_bonding_perc"] == 0.73333
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP_antibonding_perc"] == 0.26667
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["ICOBI"] == 0.08932
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["ICOOP"] == 0.0148
        assert graph_cdf_all.sg.graph.get_edge_data(0, 1)[3]["bond_label"] == "63"

    def test_graph_cdf_close_fermi(self):
        graph_cdf_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
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
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[0]["specie"] == "Cd"
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist() == [0.0, 0.0, 0.0]
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["Mulliken Charges"] == 1.57
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[0]["properties"]["env"] == "C:8"

        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[1]["specie"] == "F"
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist() == [1.37314, 1.37314, 1.37314]
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["Mulliken Charges"] == -0.79
        assert graph_cdf_close_fermi.sg.graph.nodes(data=True)[1]["properties"]["env"] == "T:4"

        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == -0.62168
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"] == 0
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP_antibonding_perc"] == 1
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == 0.08932
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == 0.0148
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "29"

        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP"] == -0.62168
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP_bonding_perc"] == 0
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP_antibonding_perc"] == 1
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOBI"] == 0.08932
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOOP"] == 0.0148
        assert graph_cdf_close_fermi.sg.graph.get_edge_data(0, 1)[3]["bond_label"] == "63"

    def test_graph_exceptions(self):
        with pytest.raises(ValueError) as err1:  # noqa: PT012, PT011
            _ = LobsterGraph(
                path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
                path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
                path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
                add_additional_data_sg=True,
                path_to_icooplist=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz",
                path_to_icobilist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
                path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
                which_bonds="cation-anin",
                start=None,
            )

            assert (
                str(err1.value) == "Only accepted values are 'all' and 'cation-anion'. "
                "Please check the input parameters of which_bonds arg"
            )

        with pytest.raises(ValueError) as err2:  # noqa: PT012, PT011
            _ = LobsterGraph(
                path_to_poscar=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
                path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
                path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
                path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
                add_additional_data_sg=True,
                path_to_icooplist=None,
                path_to_icobilist=None,
                path_to_madelung=TestDir / "test_data/NaCl_comp_range/MadelungEnergies.lobster.gz",
                which_bonds="cation-anion",
                start=None,
            )

            assert (
                str(err2.value) == "add_additional_data_sg is set to True. "
                "Please provide path_to_icooplist and path_to_icobilist"
            )
