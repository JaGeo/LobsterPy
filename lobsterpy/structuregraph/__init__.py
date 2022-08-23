# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This package provides the modules for generating graph objects using lobsterpy data
"""
from pymatgen.core.structure import Structure
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.lobster.outputs import Charge


class GraphObject:
    """
    GraphObject class for bonding data from Lobster

    Attributes:
        sg: return structure_graph object
    """

    def __init__(
        self,
        path_to_poscar: str,
        path_to_charge: str,
        path_to_icohplist: str,
        path_to_icooplist: str = None,
        path_to_icobilist: str = None,
        add_additional_data_sg=True,
    ):
        """
        This is a class to get structure graph object with bonding information

        Args:
            path_to_poscar: str that describes path to POSCAR
            path_to_charge: str that describes the path to CHARGE.lobster
            path_to_icohplist: str that describes the path to ICOHPLIST.lobster
            path_to_icooplist: str that describes the path to ICOOPLIST.lobster
            path_to_icobilist: str that describes the path to ICOBILIST.lobster
            add_additional_data_sg: (bool) will add the information from ICOOPLIST.lobster and ICOBILIST.lobster
        """

        self.path_to_poscar = path_to_poscar
        self.path_to_charge = path_to_charge
        self.path_to_icohplist = path_to_icohplist
        self.path_to_icooplist = path_to_icooplist
        self.path_to_icobilist = path_to_icobilist
        self.add_additional_data_sg = add_additional_data_sg

    def sg(self):

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

        decorated_structure = Charge(self.path_to_charge).get_structure_with_charges(
            self.path_to_poscar
        )

        # decorate=True will add the coordination environments as ell
        sg = chemenvlobster.get_bonded_structure(
            structure=decorated_structure, decorate=True, edge_properties=True
        )

        return sg

    def get_sg_node_data(self):

        graph_data = self.sg()
        print(graph_data.graph.nodes(data=True))

    def get_sg_edge_data(self):

        graph_data = self.sg()
        print(graph_data.graph.edges(data=True))
