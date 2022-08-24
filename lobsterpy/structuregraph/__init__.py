# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This package provides the modules for generating graph objects using lobsterpy data
"""

import os
from pymatgen.core.structure import Structure
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.lobster.outputs import Charge
from lobsterpy.cohp.analyze import Analysis


class GraphObject:
    """
    GraphObject class for bonding data from Lobster

    Attributes:
        sg: return structure_graph object
    """

    def __init__(
        self,
        path_to_calcdir: str,
        add_additional_data_sg=True,
    ):
        """
        This is a class to get structure graph object with bonding information from lobsterpy

        Args:
            path_to_calcdir: str that describes path to Lobster run files
            (should include POSCAR, CHARGE.lobster, ICOHPLIST.lobster, ICOOPLIST.lobster, ICOBILIST.lobster,
            MadelungEnergies.lobster)
            add_additional_data_sg: (bool) will add the information from ICOOPLIST.lobster and ICOBILIST.lobster
        """
        self.path_to_calcdir = path_to_calcdir

        if add_additional_data_sg:
            self.add_additional_data_sg = add_additional_data_sg
            self.path_to_icooplist = os.path.join(
                self.path_to_calcdir, "ICOOPLIST.lobster"
            )
            self.path_to_icobilist = os.path.join(
                self.path_to_calcdir, "ICOBILIST.lobster"
            )

        else:
            self.add_additional_data_sg = add_additional_data_sg

        self.path_to_poscar = os.path.join(self.path_to_calcdir, "POSCAR")
        self.path_to_charge = os.path.join(self.path_to_calcdir, "CHARGE.lobster")
        self.path_to_cohpcar = os.path.join(self.path_to_calcdir, "COHPCAR.lobster")
        self.path_to_icohplist = os.path.join(self.path_to_calcdir, "ICOHPLIST.lobster")
        self.path_to_madelung = os.path.join(
            self.path_to_calcdir, "MadelungEnergies.lobster"
        )
        self.sg = self.get_decorated_sg()

    def get_decorated_sg(self):
        """
        Method to generate graph object decorated with bonding data from lobsterpy
        Returns:
            structure graph object
        """
        if self.add_additional_data_sg:
            chemenvlobster = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=self.path_to_icohplist,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=0,
                filename_CHARGE=self.path_to_charge,
                add_additional_data_sg=self.add_additional_data_sg,
                filename_blist_sg1=self.path_to_icobilist,
                id_blist_sg1="ICOBI",
                filename_blist_sg2=self.path_to_icooplist,
                id_blist_sg2="ICOOP",
            )

        else:

            chemenvlobster = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=self.path_to_icohplist,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=0,
                filename_CHARGE=self.path_to_charge,
                add_additional_data_sg=self.add_additional_data_sg,
            )

        decorated_structure = Charge(self.path_to_charge).get_structure_with_charges(
            self.path_to_poscar
        )

        lobster_env = chemenvlobster.get_bonded_structure(
            structure=decorated_structure, decorate=True, edge_properties=True
        )

        analyze = Analysis(
            path_to_charge=self.path_to_charge,
            path_to_cohpcar=self.path_to_cohpcar,
            path_to_poscar=self.path_to_poscar,
            path_to_icohplist=self.path_to_icohplist,
            path_to_madelung=self.path_to_madelung,
            whichbonds="all",
        )

        cba = analyze.condensed_bonding_analysis

        for k, v in cba["sites"].items():
            for k2, v2 in lobster_env.graph.nodes.data():
                if v["ion"] == v2["specie"]:
                    v2["properties"].update({"env": v["env"]})

        return lobster_env

    def get_sg_node_data(self):
        """
        Method of extract node data of graph object
        Returns:
           will return a list consisting of node data
        """

        return print(self.sg.graph.nodes(data=True))

    def get_sg_edge_data(self):
        """
        Method of extract edge data of graph object
        Returns:
           will return a list consisting of edge data
        """

        return print(self.sg.graph.edges(data=True))
