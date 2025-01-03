# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module provides a class for generating graph objects using lobsterpy data."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from pymatgen.core.structure import Structure
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.lobster.outputs import Charge

from lobsterpy.cohp.analyze import Analysis


class LobsterGraph:
    """
    Class to generate structure graph objects with bonding data from Lobster.

    Mode of automatic bonding analysis can be `cation-anion` or `all` bonds. The strongest bond is
    determined based on the ICOHPs. The coordination environments are determined based on
    `cutoff_icohp * ICOHPs` values. If the path of ICOBILIST (ICOOPLIST) is provided, the ICOBI (ICOOP)
    values corresponding to relevant bond labels obtained from the ICOHPLIST are also added as edge properties
    to the structure graph objects. The Mulliken and Loewdin charges are added as node properties to
    the structure graph objects.

    Attributes:
        sg: return structure_graph object

    :param path_to_poscar: path to POSCAR (e.g., "POSCAR").
    :param path_to_charge: path to CHARGE.lobster (e.g., "CHARGE.lobster").
    :param path_to_cohpcar: path to COHPCAR.lobster (e.g., "COHPCAR.lobster").
    :param path_to_icohplist: path to ICOHPLIST.lobster (e.g., "ICOHPLIST.lobster").
    :param path_to_icooplist: path to ICOOPLIST.lobster (e.g., "ICOOPLIST.lobster").
    :param path_to_icobilist: path to ICOBILIST.lobster (e.g., "ICOBILIST.lobster").
    :param path_to_madelung: path to MadelungEnergies.lobster (e.g., "MadelungEnergies.lobster")
    :param cutoff_icohp: only bonds that are stronger than cutoff_icohp * strongest ICOHP will be considered.
    :param noise_cutoff: if provided hardcodes the lower limit of icohps considered.
    :param add_additional_data_sg: if True will add the information from ICOOPLIST.lobster
        and ICOBILIST.lobster based on ICOHPLIST.lobster relevant bond.
    :param which_bonds: selects which kind of bonds are analyzed. "all" is the default.
    :param start: start energy for bonding antibonding percent integration.
    """

    def __init__(
        self,
        path_to_poscar: str | Path,
        path_to_charge: str | Path,
        path_to_cohpcar: str | Path,
        path_to_icohplist: str | Path,
        path_to_madelung: str | Path,
        add_additional_data_sg: bool = True,
        path_to_icooplist: str | Path | None = None,
        path_to_icobilist: str | Path | None = None,
        which_bonds: str = "all",
        cutoff_icohp: float = 0.10,
        noise_cutoff: float = 0.1,
        start: float | None = None,
    ):
        """Initialize and return a structure graph object."""
        if add_additional_data_sg:
            self.add_additional_data_sg = add_additional_data_sg
            if path_to_icooplist is not None and path_to_icobilist is not None:
                self.path_to_icooplist = path_to_icooplist
                self.path_to_icobilist = path_to_icobilist
            else:
                raise ValueError(
                    "add_additional_data_sg is set to True.Please provide path_to_icooplist and path_to_icobilist"
                )
        else:
            self.add_additional_data_sg = add_additional_data_sg  # type: ignore

        self.path_to_poscar = path_to_poscar
        self.path_to_charge = path_to_charge
        self.path_to_cohpcar = path_to_cohpcar
        self.path_to_icohplist = path_to_icohplist
        self.path_to_madelung = path_to_madelung
        self.which_bonds = which_bonds
        self.cutoff_icohp = cutoff_icohp
        self.noise_cutoff = noise_cutoff

        if self.which_bonds == "all":
            self.additional_condition = 0
        elif self.which_bonds == "cation-anion":
            self.additional_condition = 1
        else:
            raise ValueError(
                "Only accepted values are 'all' and 'cation-anion'."
                "Please check the input parameters of which_bonds arg"
            )
        self.start = start

        self.sg = self.get_decorated_sg()

    def get_decorated_sg(self):
        """
        Create a graph object decorated with bonding data from LobsterPy.

        Returns:
            structure graph object
        """
        if self.add_additional_data_sg:
            chemenvlobster = LobsterNeighbors(
                are_coops=False,
                filename_icohp=self.path_to_icohplist,
                perc_strength_icohp=self.cutoff_icohp,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=self.additional_condition,
                filename_charge=self.path_to_charge,
                add_additional_data_sg=self.add_additional_data_sg,
                filename_blist_sg1=self.path_to_icobilist,
                id_blist_sg1="ICOBI",
                filename_blist_sg2=self.path_to_icooplist,
                id_blist_sg2="ICOOP",
                valences_from_charges=True,
                adapt_extremum_to_add_cond=True,
                noise_cutoff=self.noise_cutoff,
            )

        else:
            chemenvlobster = LobsterNeighbors(
                are_coops=False,
                filename_icohp=self.path_to_icohplist,
                perc_strength_icohp=self.cutoff_icohp,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=self.additional_condition,
                filename_charge=self.path_to_charge,
                add_additional_data_sg=self.add_additional_data_sg,
                valences_from_charges=True,
                adapt_extremum_to_add_cond=True,
                noise_cutoff=self.noise_cutoff,
            )

        # Adds Mulliken and Löwdin charges as site properties to structure object (node properties)
        decorated_structure = Charge(self.path_to_charge).get_structure_with_charges(self.path_to_poscar)

        # Create the structure graph object decorated with site and edge properties based on ICOHP/ICOBI/ICOOP data
        lobster_env = chemenvlobster.get_bonded_structure(
            structure=decorated_structure, decorate=True, edge_properties=True
        )

        # Initialize automating bonding analysis from Lobsterpy based on ICOHP
        analyze = Analysis(
            path_to_charge=self.path_to_charge,
            path_to_cohpcar=self.path_to_cohpcar,
            path_to_poscar=self.path_to_poscar,
            path_to_icohplist=self.path_to_icohplist,
            path_to_madelung=self.path_to_madelung,
            which_bonds=self.which_bonds,
            cutoff_icohp=self.cutoff_icohp,
        )

        # Store the summarized dictionary object containing bonding information
        cba = analyze.condensed_bonding_analysis

        # Iterate over sites in the dictionary
        for cba_data in cba["sites"].values():
            for _, node_data in lobster_env.graph.nodes.data():
                # Check if ions are same and add its corresponding environment in node properties
                if cba_data["ion"] == node_data["specie"]:
                    node_data["properties"].update({"env": cba_data["env"]})

        for edge_prop in lobster_env.graph.edges.data():  # Iterate over structure graph edges
            _ab, ab_p, _b, b_p = analyze._integrate_antbdstates_below_efermi(
                cohp=analyze.chemenv.completecohp.get_cohp_by_label(edge_prop[2]["bond_label"]),
                start=self.start,
            )  # Compute bonding- antibonding percentages for each bond in structure graph object
            edge_prop[2]["ICOHP_bonding_perc"] = b_p  # Store bonding percentage in edge of graph object
            edge_prop[2]["ICOHP_antibonding_perc"] = ab_p  # Store anti-bonding percentage in edge graph object

        return lobster_env
