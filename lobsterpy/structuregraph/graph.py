# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This package provides the modules for generating graph objects using lobsterpy data
"""

import os
from typing import Optional
from pymatgen.core.structure import Structure
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.lobster.outputs import Charge
from lobsterpy.cohp.analyze import Analysis

class LobsterGraph:
    """
    LobsterGraph class for bonding data from Lobster

    Attributes:
        sg: return structure_graph object
    """

    def __init__(
        self,
        path_to_poscar: str,
        path_to_charge: str,
        path_to_cohpcar: str,
        path_to_icohplist:str,
        path_to_madelung: str,
        add_additional_data_sg=True,
        path_to_icooplist: Optional[str] = None,
        path_to_icobilist: Optional[str] = None,
        which_bonds: str = 'all',
        start: str = None,
    ):
        """
        This is a class to get structure graph object with bonding information from lobsterpy

        Args:
            path_to_poscar: path to POSCAR (e.g., "POSCAR")
            path_to_charge: path to CHARGE.lobster (e.g., "CHARGE.lobster")
            path_to_cohpcar: path to COHPCAR.lobster (e.g., "COHPCAR.lobster")
            path_to_icohplist: path to ICOHPLIST.lobster (e.g., "ICOHPLIST.lobster")
            path_to_madelung: path to MadelungEnergies.lobster (e.g., "MadelungEnergies.lobster")
            add_additional_data_sg: (bool) if True will add the information from ICOOPLIST.lobster
            and ICOBILIST.lobster based on ICOHPLIST.lobster relevant bond
            whichbonds: selects which kind of bonds are analyzed. "all" is the default
            start: start energy for bonding antibonding percent integration
        """

        if add_additional_data_sg:
            self.add_additional_data_sg = add_additional_data_sg
            self.path_to_icooplist = path_to_icooplist
            self.path_to_icobilist = path_to_icobilist
        else:
            self.add_additional_data_sg = add_additional_data_sg
            self.path_to_icooplist = path_to_icooplist
            self.path_to_icobilist = path_to_icobilist

        self.path_to_poscar = path_to_poscar
        self.path_to_charge = path_to_charge
        self.path_to_cohpcar = path_to_cohpcar
        self.path_to_icohplist = path_to_icohplist
        self.path_to_madelung = path_to_madelung
        self.which_bonds = which_bonds

        if self.which_bonds == 'all':
            self.additional_condition = 0
        elif self.which_bonds == 'cation-anion':
            self.additional_condition = 1
        else:
            raise ValueError ("Only accepted values are 'all' and 'cation-anion'."
                              "Please check the input parameters of which_bonds arg")
        self.start = start

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
                additional_condition=self.additional_condition,
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
                additional_condition=self.additional_condition,
                filename_CHARGE=self.path_to_charge,
                add_additional_data_sg=self.add_additional_data_sg,
            )

        chemenvlobster.get_info_cohps_to_neighbors(path_to_COHPCAR=self.path_to_cohpcar)

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
            whichbonds=self.which_bonds,
        )

        cba = analyze.condensed_bonding_analysis

        for k, v in cba["sites"].items():
            for k2, v2 in lobster_env.graph.nodes.data():
                if v["ion"] == v2["specie"]:
                    v2["properties"].update({"env": v["env"]})

        for edge_prop in lobster_env.graph.edges.data():
            _ab, ab_p, _b, b_p = analyze._integrate_antbdstates_below_efermi(
                cohp=chemenvlobster.completecohp.get_cohp_by_label(edge_prop[2]["ICOHP_bond_key"]),
                start=self.start,
            )
            edge_prop[2]["ICOHP_bonding_perc"] = b_p
            edge_prop[2]["ICOHP_antibonding_perc"] = ab_p

        return lobster_env