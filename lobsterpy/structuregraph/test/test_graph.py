import unittest
from pathlib import Path

from lobsterpy.structuregraph.graph import LobsterGraph

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.graph_NaCl_all = LobsterGraph(
            path_to_poscar=TestDir / "TestData/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "TestData/NaCl_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "TestData/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/NaCl_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=None,
        )

        self.graph_NaCl_cation_anion = LobsterGraph(
            path_to_poscar=TestDir / "TestData/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "TestData/NaCl_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "TestData/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/NaCl_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="cation-anion",
            start=None,
        )

        self.graph_NaCl_without_add_data = LobsterGraph(
            path_to_poscar=TestDir / "TestData/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=None,
        )

        self.graph_NaCl_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "TestData/NaCl_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/NaCl_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/NaCl_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/NaCl_comp_range/MadelungEnergies.lobster.gz",
            add_additional_data_sg=False,
            which_bonds="all",
            start=-4,
        )

        self.graph_CdF_all = LobsterGraph(
            path_to_poscar=TestDir / "TestData/CdF_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/CdF_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/CdF_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "TestData/CdF_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "TestData/CdF_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/CdF_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=None,
        )

        self.graph_CdF_close_fermi = LobsterGraph(
            path_to_poscar=TestDir / "TestData/CdF_comp_range/POSCAR.gz",
            path_to_charge=TestDir / "TestData/CdF_comp_range/CHARGE.lobster.gz",
            path_to_cohpcar=TestDir / "TestData/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/CdF_comp_range/ICOHPLIST.lobster.gz",
            add_additional_data_sg=True,
            path_to_icooplist=TestDir / "TestData/CdF_comp_range/ICOOPLIST.lobster.gz",
            path_to_icobilist=TestDir / "TestData/CdF_comp_range/ICOBILIST.lobster.gz",
            path_to_madelung=TestDir
            / "TestData/CdF_comp_range/MadelungEnergies.lobster.gz",
            which_bonds="all",
            start=-2,
        )

    def test_graph_NaCl_all(self):
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["specie"], "Na"
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["coords"].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            0.78,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[0]["properties"]["env"], "O:6"
        )

        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["specie"], "Cl"
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["coords"].tolist(),
            [2.845847, 2.845847, 2.845847],
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.78,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.nodes(data=True)[1]["properties"]["env"], "O:6"
        )

        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"],
            -0.56614,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"], 0.08484
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"], 0.02826
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"], "21"
        )

        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP"], -0.56612
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOHP_bonding_perc"], 1
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOBI"], 0.08482
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["ICOOP"], 0.02824
        )
        self.assertEqual(
            self.graph_NaCl_all.sg.graph.get_edge_data(0, 1)[4]["bond_label"], "28"
        )

    def test_graph_NaCl_cation_anion(self):
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["specie"], "Na"
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0][
                "coords"
            ].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            0.78,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[0]["properties"][
                "env"
            ],
            "O:6",
        )

        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["specie"], "Cl"
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1][
                "coords"
            ].tolist(),
            [2.845847, 2.845847, 2.845847],
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.78,
        )
        self.assertNotIn(
            "env",
            self.graph_NaCl_cation_anion.sg.graph.nodes(data=True)[1]["properties"],
        )

        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOHP"],
            -0.56614,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOBI"],
            0.08484,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["ICOOP"],
            0.02826,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[0]["bond_label"],
            "21",
        )

        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOHP"],
            -0.56612,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOBI"],
            0.08482,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["ICOOP"],
            0.02824,
        )
        self.assertEqual(
            self.graph_NaCl_cation_anion.sg.graph.get_edge_data(0, 1)[4]["bond_label"],
            "28",
        )

    def test_graph_NaCl_without_add_data(self):
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["specie"],
            "Na",
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0][
                "coords"
            ].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            0.78,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[0]["properties"][
                "env"
            ],
            "O:6",
        )

        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["specie"],
            "Cl",
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1][
                "coords"
            ].tolist(),
            [2.845847, 2.845847, 2.845847],
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.78,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.nodes(data=True)[1]["properties"][
                "env"
            ],
            "O:6",
        )

        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]["ICOHP"],
            -0.56614,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertNotIn(
            "ICOBI", self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        )
        self.assertNotIn(
            "ICOOP", self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0]
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[0][
                "bond_label"
            ],
            "21",
        )

        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]["ICOHP"],
            -0.56612,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertNotIn(
            "ICOBI", self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        )
        self.assertNotIn(
            "ICOOP", self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4]
        )
        self.assertEqual(
            self.graph_NaCl_without_add_data.sg.graph.get_edge_data(0, 1)[4][
                "bond_label"
            ],
            "28",
        )

    def test_graph_NaCl_close_fermi(self):
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["specie"], "Na"
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            0.78,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[0]["properties"][
                "env"
            ],
            "O:6",
        )

        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["specie"], "Cl"
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist(),
            [2.845847, 2.845847, 2.845847],
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.78,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.nodes(data=True)[1]["properties"][
                "env"
            ],
            "O:6",
        )

        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"],
            -0.56614,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertNotIn(
            "ICOBI", self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        )
        self.assertNotIn(
            "ICOOP", self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"],
            "21",
        )

        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["ICOHP"],
            -0.56612,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_bonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4][
                "ICOHP_antibonding_perc"
            ],
            0,
        )
        self.assertNotIn(
            "ICOBI", self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        )
        self.assertNotIn(
            "ICOOP", self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]
        )
        self.assertEqual(
            self.graph_NaCl_close_fermi.sg.graph.get_edge_data(0, 1)[4]["bond_label"],
            "28",
        )

    def test_graph_CdF_all(self):
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[0]["specie"], "Cd"
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[0]["coords"].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            1.57,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[0]["properties"]["env"], "C:8"
        )

        self.assertEqual(self.graph_CdF_all.sg.graph.nodes(data=True)[1]["specie"], "F")
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[1]["coords"].tolist(),
            [1.37314, 1.37314, 1.37314],
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.79,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.nodes(data=True)[1]["properties"]["env"], "T:4"
        )

        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP"], -0.62168
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOHP_bonding_perc"],
            0.72263,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            0.27737,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOBI"], 0.08932
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["ICOOP"], 0.0148
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[0]["bond_label"], "29"
        )

        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP"], -0.62168
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOHP_bonding_perc"],
            0.72263,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3][
                "ICOHP_antibonding_perc"
            ],
            0.27737,
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOBI"], 0.08932
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["ICOOP"], 0.0148
        )
        self.assertEqual(
            self.graph_CdF_all.sg.graph.get_edge_data(0, 1)[3]["bond_label"], "63"
        )

    def test_graph_CdF_close_fermi(self):
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["specie"], "Cd"
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["coords"].tolist(),
            [0.0, 0.0, 0.0],
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["properties"][
                "Mulliken Charges"
            ],
            1.57,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[0]["properties"][
                "env"
            ],
            "C:8",
        )

        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["specie"], "F"
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["coords"].tolist(),
            [1.37314, 1.37314, 1.37314],
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["properties"][
                "Mulliken Charges"
            ],
            -0.79,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.nodes(data=True)[1]["properties"][
                "env"
            ],
            "T:4",
        )

        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOHP"],
            -0.62168,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_bonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0][
                "ICOHP_antibonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOBI"], 0.08932
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["ICOOP"], 0.0148
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[0]["bond_label"],
            "29",
        )

        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOHP"],
            -0.62168,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3][
                "ICOHP_bonding_perc"
            ],
            0,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3][
                "ICOHP_antibonding_perc"
            ],
            1,
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOBI"], 0.08932
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["ICOOP"], 0.0148
        )
        self.assertEqual(
            self.graph_CdF_close_fermi.sg.graph.get_edge_data(0, 1)[3]["bond_label"],
            "63",
        )


if __name__ == "__main__":
    unittest.main()
