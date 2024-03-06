# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module defines classes to analyze the COHPs automatically."""
from __future__ import annotations

import warnings
from collections import Counter
from pathlib import Path

import numpy as np
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import LobsterCompleteDos
from pymatgen.io.lobster import (
    Bandoverlaps,
    Charge,
    Doscar,
    Icohplist,
    Lobsterin,
    Lobsterout,
    MadelungEnergies,
)
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Analysis:
    """
    Class to automatically analyze COHP/COOP/COBI  populations from Lobster.

    :param are_cobis: bool indicating if file contains COBI/ICOBI data
    :param are_coops: bool indicating if file contains COOP/ICOOP data
    :param cutoff_icohp: Cutoff in percentage for evaluating neighbors based on ICOHP values.
        cutoff_icohp*max_icohp limits the number of considered neighbours for evaluating environments.
    :param path_to_cohpcar: path to `COHPCAR.lobster` or `COBICAR.lobster` or `COOPCAR.lobster` .
    :param path_to_charge: path to `CHARGE.lobster`.
    :param path_to_icohplist: path to `ICOHPLIST.lobster` or `ICOBILIST.lobster` or `ICOOPLIST.lobster`.
    :param path_to_poscar: path to structure (e.g., `POSCAR` or `POSCAR.lobster`)
    :param path_to_madelung: path to `MadelungEnergies.lobster`.
    :param charge_obj: pymatgen lobster.io.charge object
    :param completecohp_obj: pymatgen.electronic_structure.cohp.CompleteCohp object
    :param icohplist_obj: pymatgen lobster.io.Icohplist object
    :param madelung_obj: pymatgen lobster.io.MadelungEnergies object
    :param noise_cutoff: Sets the lower limit tolerance for ICOHPs or ICOOPs or ICOBIs considered
        in analysis.
    :param orbital_cutoff: Sets the minimum percentage for the orbital contribution considered to be
        relevant in orbital resolved analysis. (Affects only when orbital_resolved argument is set to True)
        Set it to 0 to get results of all orbitals in the detected relevant bonds. Default is to 0.05 i.e.
        only analyzes if orbital contribution is 5 % or more.
    :param orbital_resolved: bool indicating whether orbital wise analysis is performed
    :param type_charge: If no path_to_charge is given, Valences will be used. Otherwise, Mulliken charges.
            Löwdin charges cannot be selected at the moment.
    :param which_bonds: Selects kinds of bonds that are analyzed. `cation-anion` is the default.
        Alternatively, `all` bonds can also be selected. Support to other kinds of bonds will be
        added soon.
    :param summed_spins: if True, COHP `Spin.up` and `Spin.down` populations will be summed
    :param start: sets the lower limit of energy for evaluation of bonding and antibonding
        percentages below efermi. Defaults to None (i.e., all populations below efermi are included)

    Attributes:
        - condensed_bonding_analysis: dict including a summary of the most important bonding properties
        - final_dict_bonds: dict including information on ICOHPs per bond type
        - final_dict_ions: dict including information on environments of cations
        - chemenv: pymatgen.io.lobster.lobsterenv.LobsterNeighbors object
        - lse: LightStructureEnvironment from pymatgen
        - anion_types: Set of Element objects from pymatgen
        - list_equivalent_sites: list of site indices of sites that indicate which sites are equivalent
            e.g., [0 1 2 2 2] where site 0, 1, 2 indicate sites that are independent from each other
        - seq_cohps: list of cohps
        - seq_coord_ions: list of co-ordination environment strings for each cation
        - seq_equivalent_sites: seq of inequivalent sites
        - seq_ineq_ions: seq of inequivalent cations/sites in the structure
        - seq_infos_bonds (list): information on cation anion bonds (lists
            of pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo)
        - spg: space group information
        - structure: Structure object

    """

    def __init__(
        self,
        path_to_poscar: str | Path | None,
        path_to_icohplist: str | Path | None,
        path_to_cohpcar: str | Path | None,
        path_to_charge: str | Path | None = None,
        path_to_madelung: str | Path | None = None,
        icohplist_obj: Icohplist | None = None,
        completecohp_obj: CompleteCohp | None = None,
        charge_obj: Charge | None = None,
        madelung_obj: MadelungEnergies | None = None,
        are_cobis: bool = False,
        are_coops: bool = False,
        cutoff_icohp: float = 0.1,
        noise_cutoff: float = 0.1,
        orbital_cutoff: float = 0.05,
        orbital_resolved: bool = False,
        start: float | None = None,
        summed_spins: bool = True,
        type_charge: str | None = None,
        which_bonds: str = "cation-anion",
    ):
        """
        Initialize automatic bonding analysis.

        :param are_cobis: bool indicating if file contains COBI/ICOBI data
        :param are_coops: bool indicating if file contains COOP/ICOOP data
        :param cutoff_icohp: Cutoff in percentage for evaluating neighbors based on ICOHP values.
            cutoff_icohp*max_icohp limits the number of considered neighbours for evaluating environments.
        :param path_to_cohpcar: path to `COHPCAR.lobster` or `COBICAR.lobster` or `COOPCAR.lobster` .
        :param path_to_charge: path to `CHARGE.lobster`.
        :param path_to_icohplist: path to `ICOHPLIST.lobster` or `ICOBILIST.lobster` or `ICOOPLIST.lobster`.
        :param path_to_poscar: path to structure (e.g., `POSCAR` or `POSCAR.lobster`)
        :param path_to_madelung: path to `MadelungEnergies.lobster`.
        :param charge_obj: pymatgen lobster.io.charge object (Optional)
        :param completecohp_obj: pymatgen.electronic_structure.cohp.CompleteCohp object
        :param icohplist_obj: pymatgen lobster.io.Icohplist object
        :param madelung_obj: pymatgen lobster.io.MadelungEnergies object
        :param noise_cutoff: Sets the lower limit tolerance for ICOHPs or ICOOPs or ICOBIs considered
            in analysis.
        :param orbital_cutoff: Sets the minimum percentage for the orbital contribution considered to be
            relevant in orbital resolved analysis. (Affects only when orbital_resolved argument is set to True)
            Set it to 0 to get results of all orbitals in the detected relevant bonds. Default is to 0.05 i.e.
            only analyzes if orbital contribution is 5 % or more.
        :param orbital_resolved: bool indicating whether orbital wise analysis is performed
        :param type_charge: If no path_to_charge is given, Valences will be used. Otherwise, Mulliken charges.
                Löwdin charges cannot be selected at the moment.
        :param which_bonds: Selects kinds of bonds that are analyzed. `cation-anion` is the default.
            Alternatively, `all` bonds can also be selected. Support to other kinds of bonds will be
            added soon.
        :param summed_spins: if True, COHP `Spin.up` and `Spin.down` populations will be summed
        :param start: sets the lower limit of energy for evaluation of bonding and antibonding
            percentages below efermi. Defaults to None (i.e., all populations below efermi are included)

        """
        self.start = start
        self.completecohp_obj = completecohp_obj
        self.icohplist_obj = icohplist_obj
        # checks to ensure LobsterEnv inputs are not duplicated in case users provide both path and obj
        if self.completecohp_obj is not None and self.icohplist_obj is not None:
            self.path_to_poscar = None
            self.path_to_cohpcar = None
            self.path_to_icohplist = None
        else:
            self.path_to_poscar = path_to_poscar
            self.path_to_icohplist = path_to_icohplist
            self.path_to_cohpcar = path_to_cohpcar
        self.which_bonds = which_bonds
        self.cutoff_icohp = cutoff_icohp
        self.orbital_cutoff = orbital_cutoff
        self.path_to_charge = path_to_charge
        self.charge_obj = charge_obj
        self.path_to_madelung = path_to_madelung
        self.madelung_obj = madelung_obj
        self.are_cobis = are_cobis
        self.are_coops = are_coops
        self.noise_cutoff = noise_cutoff
        self.setup_env()
        self.get_information_all_bonds(summed_spins=summed_spins)
        self.orbital_resolved = orbital_resolved

        # This determines how cations and anions
        if path_to_charge is None and charge_obj is None:
            self.type_charge = "Valences"
        else:
            if type_charge is None:
                self.type_charge = "Mulliken"
            elif type_charge == "Mulliken":
                self.type_charge = "Mulliken"
            elif type_charge == "Löwdin":
                raise ValueError("Only Mulliken charges can be used here at the moment. Implementation will follow.")
            else:
                self.type_charge = "Valences"
                print("type_charge cannot be read! Please use Mulliken/Löwdin. Now, we will use valences")

        self.set_condensed_bonding_analysis()
        self.set_summary_dicts()
        self.path_to_madelung = path_to_madelung

    def setup_env(self):
        """
        Set up the light structure environments based on COHPs using this method.

        Returns:
            None

        """
        self.structure = (
            Structure.from_file(self.path_to_poscar) if self.path_to_poscar else self.completecohp_obj.structure
        )
        sga = SpacegroupAnalyzer(structure=self.structure)
        symmetry_dataset = sga.get_symmetry_dataset()
        equivalent_sites = symmetry_dataset["equivalent_atoms"]
        self.list_equivalent_sites = equivalent_sites
        self.seq_equivalent_sites = list(set(equivalent_sites))
        self.spg = symmetry_dataset["international"]

        if self.which_bonds == "cation-anion":
            try:
                self.chemenv = LobsterNeighbors(
                    filename_icohp=self.path_to_icohplist,
                    obj_icohp=self.icohplist_obj,
                    structure=self.structure,
                    additional_condition=1,
                    perc_strength_icohp=self.cutoff_icohp,
                    filename_charge=self.path_to_charge,
                    obj_charge=self.charge_obj,
                    valences=None,
                    valences_from_charges=True,
                    adapt_extremum_to_add_cond=True,
                    are_cobis=self.are_cobis,
                    are_coops=self.are_coops,
                    noise_cutoff=self.noise_cutoff,
                )
            except ValueError as err:
                if (
                    str(err) == "min() arg is an empty sequence"
                    or str(err) == "All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work"
                ):
                    raise ValueError(
                        "Consider switching to an analysis of all bonds and not only cation-anion bonds."
                        " It looks like no cations are detected."
                    )
                raise err
        elif self.which_bonds == "all":
            # raise ValueError("only cation anion bonds implemented so far")
            self.chemenv = LobsterNeighbors(
                filename_icohp=self.path_to_icohplist,
                obj_icohp=self.icohplist_obj,
                structure=self.structure,
                additional_condition=0,
                perc_strength_icohp=self.cutoff_icohp,
                filename_charge=self.path_to_charge,
                obj_charge=self.charge_obj,
                valences_from_charges=True,
                adapt_extremum_to_add_cond=True,
                are_cobis=self.are_cobis,
                are_coops=self.are_coops,
                noise_cutoff=self.noise_cutoff,
            )

        else:
            raise ValueError("only cation anion and all bonds implemented so far")

        # determine cations and anions
        try:
            if self.which_bonds == "cation-anion":
                self.lse = self.chemenv.get_light_structure_environment(only_cation_environments=True)
            elif self.which_bonds == "all":
                self.lse = self.chemenv.get_light_structure_environment(only_cation_environments=False)
        except ValueError:

            class Lse:
                """Test class when error was raised."""

                def __init__(self, chemenv, valences=None):
                    """
                    Test class when error was raised.

                    :param chemenv: LobsterNeighbors object
                    :param valences: list of valences

                    """
                    if valences is None:
                        self.coordination_environments = []
                        for coord in chemenv:
                            if len(coord) > 0:
                                self.coordination_environments.append([{"ce_symbol": str(len(coord))}])
                            else:
                                self.coordination_environments.append([{"ce_symbol": None}])
                    else:
                        self.coordination_environments = []

                        for val, coord in zip(valences, chemenv):
                            if val >= 0.0 and len(coord) > 0:
                                self.coordination_environments.append([{"ce_symbol": str(len(coord))}])
                            else:
                                self.coordination_environments.append([{"ce_symbol": None}])

            if self.which_bonds == "all":
                self.lse = Lse(self.chemenv.list_coords)
            elif self.which_bonds == "cation-anion":
                # make a new list
                self.lse = Lse(self.chemenv.list_coords, self.chemenv.valences)

    def get_information_all_bonds(self, summed_spins: bool = True):
        """
        Gather all information on the bonds within the compound with this method.

        Returns:
            None

        """
        if self.which_bonds == "cation-anion":
            # this will only analyze cation anion bonds which simplifies the analysis
            self.seq_ineq_ions = []
            self.seq_coord_ions = []
            self.seq_infos_bonds = []
            self.seq_labels_cohps = []
            self.seq_cohps = []
            # only_bonds_to

            self.anion_types = self.chemenv.anion_types
            for ice, ce in enumerate(self.lse.coordination_environments):
                # only look at inequivalent sites (use of symmetry to speed everything up!)!
                # only look at those cations that have cation-anion bonds
                if ice in self.seq_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.seq_ineq_ions.append(ice)

                    self.seq_coord_ions.append(ce[0]["ce_symbol"])
                    self.seq_infos_bonds.append(self.chemenv.get_info_icohps_to_neighbors([ice]))

                    aniontype_labels = []
                    aniontype_cohps = []

                    # go through all anions in the structure!
                    for anion in self.anion_types:
                        # get labels and summed cohp objects
                        labels, summedcohps = self.chemenv.get_info_cohps_to_neighbors(
                            path_to_cohpcar=self.path_to_cohpcar,
                            obj_cohpcar=self.completecohp_obj,
                            isites=[ice],
                            summed_spin_channels=summed_spins,
                            per_bond=False,
                            only_bonds_to=[str(anion)],
                        )

                        aniontype_labels.append(labels)
                        aniontype_cohps.append(summedcohps)

                    self.seq_labels_cohps.append(aniontype_labels)
                    self.seq_cohps.append(aniontype_cohps)

        elif self.which_bonds == "all":
            # this will only analyze all bonds

            self.seq_ineq_ions = []
            self.seq_coord_ions = []
            self.seq_infos_bonds = []
            self.seq_labels_cohps = []
            self.seq_cohps = []
            # only_bonds_to
            self.elements = self.structure.composition.elements
            # self.anion_types = self.chemenv.get_anion_types()
            for ice, ce in enumerate(self.lse.coordination_environments):
                # only look at inequivalent sites (use of symmetry to speed everything up!)!
                # only look at those cations that have cation-anion bonds
                if ice in self.seq_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.seq_ineq_ions.append(ice)
                    self.seq_coord_ions.append(ce[0]["ce_symbol"])
                    self.seq_infos_bonds.append(self.chemenv.get_info_icohps_to_neighbors([ice]))

                    type_labels = []
                    type_cohps = []

                    for element in self.elements:
                        # get labels and summed cohp objects
                        labels, summedcohps = self.chemenv.get_info_cohps_to_neighbors(
                            path_to_cohpcar=self.path_to_cohpcar,
                            obj_cohpcar=self.completecohp_obj,
                            isites=[ice],
                            onlycation_isites=False,
                            summed_spin_channels=summed_spins,
                            per_bond=False,
                            only_bonds_to=[str(element)],
                        )

                        type_labels.append(labels)
                        type_cohps.append(summedcohps)

                    self.seq_labels_cohps.append(type_labels)
                    self.seq_cohps.append(type_cohps)

    def get_site_bond_resolved_labels(self):
        """
        Return relevant bond labels for each symmetrically independent site.

        Returns:
            dict with bond labels for each site, e.g.
            {'Na1: Na-Cl': ['21', '23', '24', '27', '28', '30']}

        """
        bonds = [[] for _ in range(len(self.seq_infos_bonds))]  # type: ignore
        labels = [[] for _ in range(len(self.seq_infos_bonds))]  # type: ignore
        for inx, bond_info in enumerate(self.seq_infos_bonds):
            for ixx, val in enumerate(bond_info.atoms):
                label_srt = sorted(val.copy())
                bonds[inx].append(
                    self.structure.sites[bond_info.central_isites[0]].species_string
                    + str(bond_info.central_isites[0] + 1)
                    + ": "
                    + label_srt[0].strip("0123456789")
                    + "-"
                    + label_srt[1].strip("0123456789")
                )
                labels[inx].append(bond_info.labels[ixx])

        label_data = {}
        for indx, atom_pairs in enumerate(bonds):
            searched_atom_pairs = set(atom_pairs)
            for search_item in searched_atom_pairs:
                indices = [i for i, pair in enumerate(atom_pairs) if pair == search_item]
                filtered_bond_label_list = [labels[indx][i] for i in indices]
                label_data.update({search_item: filtered_bond_label_list})

        return label_data

    def _get_orbital_resolved_data(
        self,
        nameion: str,
        iion: int,
        labels: list[str],
        bond_resolved_labels: dict[str, list[str]],
        type_pop: str,
    ):
        """
        Retrieve orbital-wise analysis data.

        :param nameion: name of symmetrically relevant cation or anion
        :param iion: index of symmetrically relevant cation or anion
        :param labels: list of bond label names
        :param bond_resolved_labels: dict of bond labels from ICOHPLIST resolved for each bond
        :param type_pop: population type analyzed. e.g. COHP or COOP or COBI

        Returns:
            dict consisting of relevant orbitals (contribution > 5 % to overall ICOHP or ICOBI or ICOOP),
            bonding and antibonding percentages with bond label names as keys.
        """
        orb_resolved_bond_info = {}
        for label in labels:
            if label is not None:
                bond_resolved_label_key = nameion + str(iion + 1) + ":" + label.split("x")[-1]
                bond_labels = bond_resolved_labels[bond_resolved_label_key]
                available_orbitals = list(self.chemenv.completecohp.orb_res_cohp[bond_labels[0]].keys())
                # initialize empty list to store orb paris for bonding,
                # antibonding integrals and  percentages
                bndg_orb_pair_list = []
                bndg_orb_integral_list = []
                bndg_orb_perc_list = []
                bndg_orb_icohp_list = []
                antibndg_orb_integral_list = []
                antibndg_orb_perc_list = []
                antibndg_orb_pair_list = []
                antibndg_orb_icohp_list = []

                # get total summed cohps using label list
                cohp_summed = self.chemenv.completecohp.get_summed_cohp_by_label_list(label_list=bond_labels)
                if type_pop.lower() == "cohp":
                    (
                        antibndg_tot,
                        per_anti_tot,
                        bndg_tot,
                        per_bndg_tot,
                    ) = self._integrate_antbdstates_below_efermi(cohp=cohp_summed, start=self.start)
                else:
                    (
                        bndg_tot,
                        per_bndg_tot,
                        antibndg_tot,
                        per_anti_tot,
                    ) = self._integrate_antbdstates_below_efermi(cohp=cohp_summed, start=self.start)

                orb_bonding_dict_data = {}  # type: ignore
                # For each orbital collect the contributions of summed bonding
                # and antibonding interactions separately
                for orb in available_orbitals:
                    cohp_summed_orb = self.chemenv.completecohp.get_summed_cohp_by_label_and_orbital_list(
                        label_list=bond_labels, orbital_list=[orb] * len(bond_labels)
                    )

                    if type_pop.lower() == "cohp":
                        (
                            antibndg_orb,
                            per_anti_orb,
                            bndg_orb,
                            per_bndg_orb,
                        ) = self._integrate_antbdstates_below_efermi(cohp=cohp_summed_orb, start=self.start)
                    else:
                        (
                            bndg_orb,
                            per_bndg_orb,
                            antibndg_orb,
                            per_anti_orb,
                        ) = self._integrate_antbdstates_below_efermi(cohp=cohp_summed_orb, start=self.start)

                    # replace nan values with zero (tackle numerical integration issues)
                    bndg_orb = bndg_orb if not np.isnan(bndg_orb) else 0
                    per_bndg_orb = per_bndg_orb if not np.isnan(per_bndg_orb) else 0
                    bndg_tot = bndg_tot if not np.isnan(bndg_tot) else 0
                    per_bndg_tot = per_bndg_tot if not np.isnan(per_bndg_tot) else 0
                    # skip collecting orb contributions if no summed bonding contribution exists
                    if bndg_tot > 0:
                        orb_icohps_bndg = []
                        for bond_label in bond_labels:
                            orb_icohp_bn = self.chemenv.Icohpcollection.get_icohp_by_label(
                                label=bond_label, orbitals=orb
                            )
                            orb_icohps_bndg.append(orb_icohp_bn)
                        bndg_orb_pair_list.append(orb)
                        bndg_orb_icohp_list.append(orb_icohps_bndg)
                        bndg_orb_integral_list.append(bndg_orb)
                        bndg_orb_perc_list.append(per_bndg_orb)

                    # replace nan values with zero (tackle numerical integration issues)
                    antibndg_orb = antibndg_orb if not np.isnan(antibndg_orb) else 0
                    per_anti_orb = per_anti_orb if not np.isnan(per_anti_orb) else 0
                    antibndg_tot = antibndg_tot if not np.isnan(antibndg_tot) else 0
                    per_anti_tot = per_anti_tot if not np.isnan(per_anti_tot) else 0
                    # skip collecting orb contributions if no summed antibonding contribution exists
                    if antibndg_tot > 0:
                        orb_icohps_anti = []
                        for bond_label in bond_labels:
                            orb_icohp_an = self.chemenv.Icohpcollection.get_icohp_by_label(
                                label=bond_label, orbitals=orb
                            )
                            orb_icohps_anti.append(orb_icohp_an)

                        antibndg_orb_pair_list.append(orb)
                        antibndg_orb_icohp_list.append(orb_icohps_anti)
                        antibndg_orb_integral_list.append(antibndg_orb)
                        antibndg_orb_perc_list.append(per_anti_orb)

                # Populate the dictionary with relevant orbitals for bonding interactions
                for inx, bndg_orb_pair in enumerate(bndg_orb_pair_list):
                    bndg_contri_perc = round(bndg_orb_integral_list[inx] / sum(bndg_orb_integral_list), 2)
                    # filter out very small bonding interactions (<self.orbital_cutoff)
                    if bndg_contri_perc > self.orbital_cutoff:
                        if bndg_orb_pair in orb_bonding_dict_data:
                            orb_bonding_dict_data[bndg_orb_pair].update(
                                {
                                    "orb_contribution_perc_bonding": bndg_contri_perc,
                                    "bonding": {
                                        "integral": bndg_orb_integral_list[inx],
                                        "perc": bndg_orb_perc_list[inx],
                                    },
                                }
                            )
                        else:
                            orb_bonding_dict_data[bndg_orb_pair] = {
                                f"I{type_pop}_mean": round(np.mean(bndg_orb_icohp_list[inx]), 4),
                                f"I{type_pop}_sum": round(np.sum(bndg_orb_icohp_list[inx]), 4),
                                "orb_contribution_perc_bonding": round(
                                    bndg_orb_integral_list[inx] / sum(bndg_orb_integral_list),
                                    2,
                                ),
                                "bonding": {
                                    "integral": bndg_orb_integral_list[inx],
                                    "perc": bndg_orb_perc_list[inx],
                                },
                            }

                # Populate the dictionary with relevant orbitals for antibonding interactions
                for inx, antibndg_orb_pair in enumerate(antibndg_orb_pair_list):
                    antibndg_contri_perc = round(
                        antibndg_orb_integral_list[inx] / sum(antibndg_orb_integral_list),
                        2,
                    )
                    # filter out very small antibonding interactions (<self.orbital_cutoff)
                    if antibndg_contri_perc > self.orbital_cutoff:
                        if antibndg_orb_pair in orb_bonding_dict_data:
                            orb_bonding_dict_data[antibndg_orb_pair].update(
                                {
                                    "orb_contribution_perc_antibonding": round(
                                        antibndg_orb_integral_list[inx] / sum(antibndg_orb_integral_list),
                                        2,
                                    ),
                                    "antibonding": {
                                        "integral": antibndg_orb_integral_list[inx],
                                        "perc": antibndg_orb_perc_list[inx],
                                    },
                                }
                            )
                        else:
                            orb_bonding_dict_data[antibndg_orb_pair] = {
                                f"I{type_pop}_mean": round(np.mean(antibndg_orb_icohp_list[inx]), 4),
                                f"I{type_pop}_sum": round(np.sum(antibndg_orb_icohp_list[inx]), 4),
                                "orb_contribution_perc_antibonding": round(
                                    antibndg_orb_integral_list[inx] / sum(antibndg_orb_integral_list),
                                    2,
                                ),
                                "antibonding": {
                                    "integral": antibndg_orb_integral_list[inx],
                                    "perc": antibndg_orb_perc_list[inx],
                                },
                            }

                orb_bonding_dict_data["relevant_bonds"] = bond_labels  # type: ignore

                orb_resolved_bond_info[bond_resolved_label_key] = orb_bonding_dict_data

        return orb_resolved_bond_info

    def _get_bond_resolved_data_stats(self, orb_resolved_bond_data: dict):
        """
        Retrieve the maximum bonding and anti-bonding orbital contributions.

        :param orb_resolved_bond_data: A dictionary with orbital names as keys and corresponding bonding data

        Returns:
            dict with orbital data stats the site for relevant orbitals, e.g.
            {'orbital_summary_stats': {'max_bonding_contribution': {'2s-3s': 0.68},
            'max_antibonding_contribution': {'2s-2pz': 0.36}}}

        """
        # get max orbital bonding and contribution for the site
        orb_pairs_bndg = []
        orb_pairs_antibndg = []
        orb_contri_bndg = []
        orb_contri_antibndg = []
        orbital_summary_stats = {"orbital_summary_stats": {}}  # type: ignore
        if orb_resolved_bond_data:
            for orb_pair, data in orb_resolved_bond_data.items():
                if "orb_contribution_perc_bonding" in data:
                    orb_pairs_bndg.append(orb_pair)
                    orb_contri_bndg.append(data["orb_contribution_perc_bonding"])

                if "orb_contribution_perc_antibonding" in data:
                    orb_pairs_antibndg.append(orb_pair)
                    orb_contri_antibndg.append(data["orb_contribution_perc_antibonding"])

            if orb_contri_bndg:
                max_orb_contri_bndg = max(orb_contri_bndg)
                max_orb_contri_bndg_inxs = [
                    inx for inx, orb_contri in enumerate(orb_contri_bndg) if orb_contri == max_orb_contri_bndg
                ]
                max_orb_contri_bndg_dict = {}
                for inx in max_orb_contri_bndg_inxs:
                    max_orb_contri_bndg_dict[orb_pairs_bndg[inx]] = orb_contri_bndg[inx]
                orbital_summary_stats["orbital_summary_stats"]["max_bonding_contribution"] = max_orb_contri_bndg_dict
            if orb_contri_antibndg:
                max_orb_contri_antibndg = max(orb_contri_antibndg)
                max_antibndg_contri_inxs = [
                    inx
                    for inx, orb_anti_per in enumerate(orb_contri_antibndg)
                    if orb_anti_per == max_orb_contri_antibndg
                ]
                max_antibndg_contri_dict = {}
                for inx in max_antibndg_contri_inxs:
                    max_antibndg_contri_dict[orb_pairs_antibndg[inx]] = orb_contri_antibndg[inx]
                orbital_summary_stats["orbital_summary_stats"][
                    "max_antibonding_contribution"
                ] = max_antibndg_contri_dict

        return orbital_summary_stats

    def get_site_orbital_resolved_labels(self):
        """
        Return relevant orbitals and bond labels for each symmetrically independent site.

        Returns:
            dict with bond labels for each site for relevant orbitals, e.g.
            {'Na1: Na-Cl': {'3s-3s': ['21', '23', '24', '27', '28', '30']}

        """
        site_bond_labels = self.get_site_bond_resolved_labels()
        orb_plot_data = {atom_pair: {} for atom_pair in site_bond_labels}
        if self.orbital_resolved:
            for site_index, cba_data in self.condensed_bonding_analysis["sites"].items():
                for atom in cba_data["bonds"]:
                    for orb_pair in cba_data["bonds"][atom]["orbital_data"]:
                        if orb_pair not in ("orbital_summary_stats", "relevant_bonds"):
                            atom_pair = [cba_data["ion"], atom]
                            atom_pair.sort()
                            key = (
                                self.structure.sites[site_index].species_string
                                + str(site_index + 1)
                                + ": "
                                + "-".join(atom_pair)
                            )
                            label_list = site_bond_labels[key]
                            orb_plot_data[key].update({orb_pair: label_list})
        else:
            print("Please set orbital_resolved to True when instantiating Analysis object, to get this data")

        return orb_plot_data

    @staticmethod
    def _get_strenghts_for_each_bond(pairs: list[list[str]], strengths: list[float], nameion: str | None = None):
        """
        Return a dictionary of bond strengths.

        :param pairs: list of list including labels for the atoms, e.g., [['O3', 'Cu1'], ['O3', 'Cu1']]
        :param strengths: list that gives the icohp strengths as a float, [-1.86287, -1.86288]
        :param nameion: string including the name of the cation in the list, e.g Cu1

        Returns:
            dict including inormation on icohps for each bond type, e.g.
            {'Yb-Sb': [-1.59769, -2.14723, -1.7925, -1.60773, -1.80149, -2.14335]}


        """
        dict_strenghts = {}  # type: ignore

        for pair, strength in zip(pairs, strengths):
            if nameion is not None:
                new = [
                    LobsterNeighbors._split_string(pair[0])[0],
                    LobsterNeighbors._split_string(pair[1])[0],
                ]
                new = Analysis._sort_name(new, nameion)
                string_here = new[0] + "-" + new[1]
            else:
                new = sorted(
                    [
                        LobsterNeighbors._split_string(pair[0])[0],
                        LobsterNeighbors._split_string(pair[1])[0],
                    ]
                )
                string_here = new[0] + "-" + new[1]

            if string_here not in dict_strenghts:
                dict_strenghts[string_here] = []
            dict_strenghts[string_here].append(strength)
        return dict_strenghts

    @staticmethod
    def _sort_name(pair: list[str], nameion: str | None = None):
        """
        Place the cation first in a list of name strings.

        :param pair: ["O","Cu"]
        :param nameion: "Cu"

        Returns:
            will return list of str, e.g. ["Cu", "O"]

        """
        if nameion is not None:
            new = []
            if pair[0] == nameion:
                new.append(pair[0])
                new.append(pair[1])

            elif pair[1] == nameion:
                new.append(pair[1])
                new.append(pair[0])

        return new

    @staticmethod
    def _sort_orbital_atom_pair(
        atom_pair: list[str],
        label: str,
        complete_cohp: CompleteCohp,
        orb_pair: str,
    ):
        """
        Place the cation first in a list of name strings and add the associated orbital name alongside the atom name.

        :param atom_pair: list of atom pair with cation first eg., ["Cl","Na"]
        :param label: LOBSTER relevant bond label eg ., "3"
        :param complete_cohp: pymatgen CompleteCohp object
        :param orb_pair: relevant orbital pair eg., "2px-3s"

        Returns:
            will return list of str, e.g. ["Na(2px)", "Cl(3s)"]

        """
        orb_atom = {}  # type: ignore
        orb_pair_list = orb_pair.split("-")
        # get orbital associated to the atom and store in a dict
        for _inx, (site, site_orb) in enumerate(zip(complete_cohp.bonds[label]["sites"], orb_pair_list)):
            if site.species_string in orb_atom:  # check necessary for bonds between same atoms
                orb_atom[site.species_string].append(site_orb)
            else:
                orb_atom[site.species_string] = [site_orb]

        orb_atom_list = []
        # add orbital name next to atom_pair
        for inx, atom in enumerate(atom_pair):
            # check to ensure getting 2nd orbital if bond is between same atomic species
            if inx == 1 and len(orb_atom.get(atom)) > 1:  # type: ignore
                atom_with_orb_name = f"{atom}({orb_atom.get(atom)[1]})"  # type: ignore
            else:
                atom_with_orb_name = f"{atom}({orb_atom.get(atom)[0]})"  # type: ignore
            orb_atom_list.append(atom_with_orb_name)

        return orb_atom_list

    def _get_antibdg_states(self, cohps, labels: list[str], nameion: str | None = None, limit=0.01):
        """
        Return a dictionary containing information on anti-bonding states.

        e.g., similar to: {'Cu-O': True, 'Cu-F': True}

        :param cohps: list of pymatgen.electronic_structure.cohp.Cohp objects
        :param labels: ['2 x Cu-O', '4 x Cu-F']
        :param nameion: string of the cation name, e.g. "Cu"
        :param limit: limit to detect antibonding states

        Returns:
            dict including in formation on whether antibonding interactions exist,
            e.g., {'Cu-O': True, 'Cu-F': True}


        """
        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            if label is not None:
                if nameion is not None:
                    new = label.split(" ")[2].split("-")
                    sorted_new = self._sort_name(new, nameion)
                    new_label = sorted_new[0] + "-" + sorted_new[1]
                else:
                    new = label.split(" ")[2].split("-")
                    sorted_new = sorted(new.copy())
                    new_label = sorted_new[0] + "-" + sorted_new[1]

                antbd = cohp.has_antibnd_states_below_efermi(limit=limit)
                if Spin.down in antbd:
                    dict_antibd[new_label] = antbd[Spin.up] or antbd[Spin.down]
                else:
                    dict_antibd[new_label] = antbd[Spin.up]

        return dict_antibd

    def _integrate_antbdstates_below_efermi_for_set_cohps(self, labels: list[str], cohps, nameion: str):
        """
        Return a dictionary containing information on antibonding states.

        .. warning:: NEEDS MORE TESTS

        It is important to note that only the energy range that has been computed can be considered
        (i.e., this might not be all)

        e.g. output: {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995},
        'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}

        :param cohps: list of pymatgen.electronic_structure.cohp.Cohp objects
        :param labels: ['2 x Cu-O', '4 x Cu-F']
        :param nameion: string of the cation name, e.g. "Cu"

        Returns:
            dict including in formation on whether antibonding interactions exist,
            e.g., {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995},
            'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}}
        """
        dict_bd_antibd = {}
        for label, cohp in zip(labels, cohps):
            if label is not None:
                new = label.split(" ")[2].split("-")
                sorted_new = self._sort_name(new, nameion)
                new_label = sorted_new[0] + "-" + sorted_new[1]
                if not self.are_cobis and not self.are_coops:
                    (
                        integral,
                        perc,
                        integral2,
                        perc2,
                    ) = self._integrate_antbdstates_below_efermi(cohp, start=self.start)
                else:
                    (
                        integral2,
                        perc2,
                        integral,
                        perc,
                    ) = self._integrate_antbdstates_below_efermi(cohp, start=self.start)

                if integral == 0 and integral2 != 0.0:
                    dict_bd_antibd[new_label] = {
                        "bonding": {"integral": integral2, "perc": perc2},
                        "antibonding": {"integral": integral, "perc": 0.0},
                    }
                elif integral2 == 0.0 and integral != 0.0:
                    dict_bd_antibd[new_label] = {
                        "bonding": {"integral": integral2, "perc": 0.0},
                        "antibonding": {"integral": integral, "perc": perc},
                    }
                elif integral == 0.0 and integral2 == 0.0:
                    dict_bd_antibd[new_label] = {
                        "bonding": {"integral": integral2, "perc": 0.0},
                        "antibonding": {"integral": integral, "perc": 0.0},
                    }
                else:
                    dict_bd_antibd[new_label] = {
                        "bonding": {"integral": integral2, "perc": perc2},
                        "antibonding": {"integral": integral, "perc": perc},
                    }

        return dict_bd_antibd

    def _integrate_antbdstates_below_efermi(self, cohp, start: float | None):
        """
        Integrate the cohp data to compute bonding and anti-bonding contribution below efermi.

        .. warning:: NEEDS MORE TESTS

        This integrates the whole COHP curve that has been computed.
        The energy range is very important.
        At present the energy range considered is dependent on COHPstartEnergy
        set during lobster runs. The bonding / antibonding integral values are sensitive to this parameter.
        If COHPstartEnergy value does not cover entire range of VASP calculations then
        absolute value of ICOHP_sum might not be equivalent to (bonding- antibonding) integral values.

        :param cohp: cohp object
        :param start: integration start energy in eV , eg start = -15

        Returns:
            absolute value of antibonding, percentage value of antibonding,
            absolute value of bonding and percentage value of bonding interactions
        """
        warnings.warn(
            "The bonding, antibonding integral/percent values are numerical estimate."
            " These values are sensitive to COHPstartEnergy parameter."
            " If COHPstartEnergy value does not cover entire range of VASP calculations then"
            " absolute value of ICOHP_sum might not be equivalent to (bonding- antibonding) integral values."
        )

        from scipy.integrate import trapezoid

        def integrate_positive(y, x):
            """
            Integrate only bonding interactions of COHPs.

            :param y: COHP values
            :param x: Energy values

            Returns:
                integrated value of bonding interactions
            """
            y = np.asanyarray(y)
            x = np.asanyarray(x)

            bonding = trapezoid(y, x)

            return np.round(bonding, 2)

        def integrate_negative(y, x):
            """
            Integrate only anti-bonding interactions of COHPs.

            :param y: COHP values
            :param x: Energy values

            Returns:
                integrated value of anti-bonding interactions
            """
            y = np.asanyarray(y)
            x = np.asanyarray(x)
            antibonding = trapezoid(y, x)

            return np.round(antibonding, 2)

        # will integrate spin.up and spin.down only below efermi
        energies_corrected = cohp.energies - cohp.efermi
        summedcohp = cohp.cohp[Spin.up] + cohp.cohp[Spin.down] if Spin.down in cohp.cohp else cohp.cohp[Spin.up]

        cohp_bf = []
        en_bf = []

        for i, en in enumerate(energies_corrected):
            if (start is None) and en <= 0:
                en_bf.append(en)
                cohp_bf.append(-1 * summedcohp[i])
            if (start is not None) and 0 >= en >= start:
                en_bf.append(en)
                cohp_bf.append(-1 * summedcohp[i])

        # Separate the bonding and antibonding COHP values in separate lists
        pos = []
        en_pos = []
        neg = []
        en_neg = []

        for i, scohp in enumerate(cohp_bf):
            if scohp >= 0:
                pos.append(scohp)
                en_pos.append(energies_corrected[i])

        for i, scohp in enumerate(cohp_bf):
            if scohp <= 0:
                neg.append(-1 * scohp)
                en_neg.append(energies_corrected[i])

        antibonding = integrate_negative(y=neg, x=en_neg)

        bonding = integrate_positive(y=pos, x=en_pos)

        return (
            antibonding,
            np.round(abs(antibonding) / (abs(bonding) + abs(antibonding)), 5),
            bonding,
            np.round(abs(bonding) / (abs(bonding) + abs(antibonding)), 5),
        )

    def _get_pop_type(self):
        """
        Return the type of the input population file.

        Returns:
            A String of analysed population can be COOP/COBI/COHP
        """
        if self.are_cobis:
            type_pop = "COBI"
        elif self.are_coops:
            type_pop = "COOP"
        else:
            type_pop = "COHP"

        return type_pop

    @staticmethod
    def _get_bond_dict(
        bond_strength_dict: dict,
        small_antbd_dict: dict,
        nameion: str | None = None,
        large_antbd_dict: dict | None = None,
        type_pop: str | None = None,
    ):
        """
        Return a bond_dict that contains information for each site.

        :param bond_strength_dict: dict with bond names as key and lists of bond strengths as items
        :param small_antbd_dict: dict including if there are antibonding interactions, {'Yb-Sb': False}
        :param nameion: name of the cation, e.g. Yb
        :param large_antbd_dict: will be implemented later
        :param type_pop: population type analyzed. eg. COHP

        Returns:
            Eg., if type_pop == 'COHP', will return
            dict including information on the anion (as label) and the ICOHPs in the item of the dict
            ICOHP_mean refers to the mean ICOHP in eV
            ICOHP_sum refers to the sum of the ICOHPs in eV
            has_antibdg_states_below_Efermi is True if there are antibonding interactions below Efermi
            "number_of_bonds" will count the numbers of bonds to the cation

        Example:
            {'Sb': {'ICOHP_mean': '-1.85', 'ICOHP_sum': '-11.09',
            'has_antibdg_states_below_Efermi': False, 'number_of_bonds': 6}}

        """
        bond_dict = {}

        for key, item in bond_strength_dict.items():
            if nameion is not None:
                a = key.split("-")[0]
                b = key.split("-")[1]
                if a == nameion:
                    key_here = b
                elif b == nameion:
                    key_here = a

            if large_antbd_dict is None:
                bond_dict[key_here] = {
                    f"I{type_pop}_mean": str(round(np.mean(item), 2)),
                    f"I{type_pop}_sum": str(round(np.sum(item), 2)),
                    "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                    "number_of_bonds": len(item),
                }
            else:
                bond_dict[key_here] = {
                    f"I{type_pop}_mean": str(round(np.mean(item), 2)),
                    f"I{type_pop}_sum": str(round(np.sum(item), 2)),
                    "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                    "number_of_bonds": len(item),
                    "perc_antibdg_states_below_Efermi": large_antbd_dict[key],
                }

        return bond_dict

    def set_condensed_bonding_analysis(self):
        """
        Condense the bonding analysis into a summary dictionary.

        Returns:
            None

        """
        self.condensed_bonding_analysis = {}
        # which icohps are considered
        if self.which_bonds == "cation-anion":
            limit_icohps = self.chemenv._get_limit_from_extremum(
                self.chemenv.Icohpcollection,
                self.cutoff_icohp,
                adapt_extremum_to_add_cond=True,
                additional_condition=1,
            )
        elif self.which_bonds == "all":
            limit_icohps = self.chemenv._get_limit_from_extremum(
                self.chemenv.Icohpcollection,
                self.cutoff_icohp,
                adapt_extremum_to_add_cond=True,
                additional_condition=0,
            )
            # formula of the compound
        formula = str(self.structure.composition.reduced_formula)
        # set population type
        type_pop = self._get_pop_type()
        # how many inequivalent cations are in the structure
        if self.which_bonds == "cation-anion":
            number_considered_ions = len(self.seq_ineq_ions)
        elif self.which_bonds == "all":
            number_considered_ions = len(self.seq_ineq_ions)

        # what was the maximum bond lengths that was considered
        max_bond_lengths = max(self.chemenv.Icohpcollection._list_length)

        # what are the charges for the cations in the structure
        charge_list = self.chemenv.valences

        # dictionary including bonding information for each site
        site_dict = {}
        if self.which_bonds == "cation-anion":
            for ication, ce, cation_anion_infos, labels, cohps in zip(
                self.seq_ineq_ions,
                self.seq_coord_ions,
                self.seq_infos_bonds,
                self.seq_labels_cohps,
                self.seq_cohps,
            ):
                namecation = str(self.structure[ication].specie)

                # This will compute the mean strengths of ICOHPs
                mean_icohps = self._get_strenghts_for_each_bond(
                    pairs=cation_anion_infos[4],
                    strengths=cation_anion_infos[1],
                    nameion=namecation,
                )
                # pairs, strengths, nameion
                # will collect if there are antibonding states present
                antbdg = self._get_antibdg_states(cohps, labels, namecation)
                dict_antibonding = self._integrate_antbdstates_below_efermi_for_set_cohps(
                    labels, cohps, nameion=namecation
                )
                bond_dict = self._get_bond_dict(mean_icohps, antbdg, namecation, type_pop=type_pop)
                bond_resolved_labels = self.get_site_bond_resolved_labels()

                for cation_name, icohp_data in bond_dict.items():
                    for atom_pair, bonding_data in dict_antibonding.items():
                        if namecation == atom_pair.split("-")[0] and cation_name == atom_pair.split("-")[1]:
                            icohp_data["bonding"] = bonding_data["bonding"]
                            icohp_data["antibonding"] = bonding_data["antibonding"]
                            if self.orbital_resolved:
                                # get orb resolved data to be added
                                orb_resolved_bond_info = self._get_orbital_resolved_data(
                                    nameion=namecation,
                                    iion=ication,
                                    labels=labels,
                                    bond_resolved_labels=bond_resolved_labels,
                                    type_pop=type_pop,
                                )
                                # match the dict key in bond_dict and get corresponding orbital data
                                for ion_atom_pair_orb in orb_resolved_bond_info:
                                    orb_data_atom_pair = ion_atom_pair_orb.split(": ")[-1]
                                    atom_pair_here = atom_pair.split("-")
                                    atom_pair_here.sort()
                                    if (
                                        orb_data_atom_pair == "-".join(atom_pair_here)
                                        and (namecation + str(ication + 1) + ":") in ion_atom_pair_orb
                                    ):
                                        icohp_data["orbital_data"] = orb_resolved_bond_info[ion_atom_pair_orb]

                                        orb_data_stats = self._get_bond_resolved_data_stats(
                                            orb_resolved_bond_data=orb_resolved_bond_info[ion_atom_pair_orb],
                                        )

                                        icohp_data["orbital_data"].update(orb_data_stats)

                site_dict[ication] = {
                    "env": ce,
                    "bonds": bond_dict,
                    "ion": namecation,
                    "charge": charge_list[ication],
                    "relevant_bonds": cation_anion_infos[3],
                }
        elif self.which_bonds == "all":
            for iion, ce, bond_infos, labels, cohps in zip(
                self.seq_ineq_ions,
                self.seq_coord_ions,
                self.seq_infos_bonds,
                self.seq_labels_cohps,
                self.seq_cohps,
            ):
                nameion = str(self.structure[iion].specie)

                # This will compute the mean strengths of ICOHPs
                mean_icohps = self._get_strenghts_for_each_bond(
                    pairs=bond_infos[4], strengths=bond_infos[1], nameion=None
                )
                # pairs, strengths, nameion
                # will collect if there are antibonding states present
                antbdg = self._get_antibdg_states(cohps, labels, nameion=None)

                dict_antibonding = self._integrate_antbdstates_below_efermi_for_set_cohps(labels, cohps, nameion)

                bond_dict = self._get_bond_dict(mean_icohps, antbdg, nameion=nameion, type_pop=type_pop)
                bond_resolved_labels = self.get_site_bond_resolved_labels()

                for cation_name, icohp_data in bond_dict.items():
                    for atom_pair, bonding_data in dict_antibonding.items():
                        if nameion == atom_pair.split("-")[0] and cation_name == atom_pair.split("-")[1]:
                            icohp_data["bonding"] = bonding_data["bonding"]
                            icohp_data["antibonding"] = bonding_data["antibonding"]
                            if self.orbital_resolved:
                                # get orb resolved data to be added
                                orb_resolved_bond_info = self._get_orbital_resolved_data(
                                    nameion=nameion,
                                    iion=iion,
                                    labels=labels,
                                    bond_resolved_labels=bond_resolved_labels,
                                    type_pop=type_pop,
                                )
                                # match the dict key in bond_dict and get corresponding orbital data
                                for ion_atom_pair_orb in orb_resolved_bond_info:
                                    orb_data_atom_pair = ion_atom_pair_orb.split(": ")[-1]
                                    atom_pair_here = atom_pair.split("-")
                                    atom_pair_here.sort()
                                    if (
                                        orb_data_atom_pair == "-".join(atom_pair_here)
                                        and (nameion + str(iion + 1) + ":") in ion_atom_pair_orb
                                    ):
                                        icohp_data["orbital_data"] = orb_resolved_bond_info[ion_atom_pair_orb]

                                        orb_data_stats = self._get_bond_resolved_data_stats(
                                            orb_resolved_bond_data=orb_resolved_bond_info[ion_atom_pair_orb],
                                        )

                                        icohp_data["orbital_data"].update(orb_data_stats)

                site_dict[iion] = {
                    "env": ce,
                    "bonds": bond_dict,
                    "ion": nameion,
                    "charge": charge_list[iion],
                    "relevant_bonds": bond_infos[3],
                }

        if self.path_to_madelung is None and self.madelung_obj is None:
            if self.which_bonds == "cation-anion":
                # This sets the dictionary including the most important information on the compound
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                }
            elif self.which_bonds == "all":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                }
        else:
            madelung = MadelungEnergies(self.path_to_madelung) if self.path_to_madelung else self.madelung_obj
            if self.type_charge == "Mulliken":
                madelung_energy = madelung.madelungenergies_mulliken
            elif self.type_charge == "Löwdin":
                madelung_energy = madelung.madelungenergies_loewdin
            else:
                madelung_energy = None
            # This sets the dictionary including the most important information on the compound
            if self.which_bonds == "cation-anion":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                    "madelung_energy": madelung_energy,
                }
            elif self.which_bonds == "all":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                    "madelung_energy": madelung_energy,
                }

    def set_summary_dicts(self):
        """
        Set summary dict that can be used for correlations.

        bond_dict that includes information on each bond

        "has_antbd" tells if there are antbonding states
        "ICOHP_mean" shows the mean of all ICOHPs in EV

        {'Yb-Sb': { 'has_antbdg': False, 'ICOHP_mean': -1.7448},
        'Mn-Sb': { 'has_antbdg': True, 'ICOHP_mean': -1.525}}

        a cation dict that includes all different coordination environments and counts for them
        {'Na': {'T:4': 4, 'A:2': 4}, 'Si': {'T:6': 4, 'PP:6': 4}}

        Returns:
            None

        """
        relevant_ion_ids = [isite for isite in self.list_equivalent_sites if isite in self.seq_ineq_ions]
        # set population type
        type_pop = self._get_pop_type()

        final_dict_bonds = {}
        for key in relevant_ion_ids:
            item = self.condensed_bonding_analysis["sites"][key]
            for type, properties in item["bonds"].items():
                label_list = [item["ion"], str(type)]
                new_label = sorted(label_list.copy())
                label = str(new_label[0]) + "-" + str(new_label[1])

                if label not in final_dict_bonds:
                    final_dict_bonds[label] = {
                        "number_of_bonds": int(properties["number_of_bonds"]),
                        f"I{type_pop}_sum": float(properties[f"I{type_pop}_sum"]),
                        "has_antbdg": properties["has_antibdg_states_below_Efermi"],
                    }
                else:
                    final_dict_bonds[label]["number_of_bonds"] += int(properties["number_of_bonds"])
                    final_dict_bonds[label][f"I{type_pop}_sum"] += float(properties[f"I{type_pop}_sum"])
                    final_dict_bonds[label]["has_antbdg"] = (
                        final_dict_bonds[label]["has_antbdg"] or properties["has_antibdg_states_below_Efermi"]
                    )
        self.final_dict_bonds = {}
        for key, item in final_dict_bonds.items():
            self.final_dict_bonds[key] = {}
            self.final_dict_bonds[key][f"I{type_pop}_mean"] = item[f"I{type_pop}_sum"] / (item["number_of_bonds"])
            self.final_dict_bonds[key]["has_antbdg"] = item["has_antbdg"]

        # rework, add all environments!
        final_dict_ions = {}
        for key in relevant_ion_ids:
            if self.condensed_bonding_analysis["sites"][key]["ion"] not in final_dict_ions:
                final_dict_ions[self.condensed_bonding_analysis["sites"][key]["ion"]] = [
                    self.condensed_bonding_analysis["sites"][key]["env"]
                ]
            else:
                final_dict_ions[self.condensed_bonding_analysis["sites"][key]["ion"]].append(
                    self.condensed_bonding_analysis["sites"][key]["env"]
                )

        self.final_dict_ions = {}
        for key, item in final_dict_ions.items():
            self.final_dict_ions[key] = dict(Counter(item))

    @staticmethod
    def get_lobster_calc_quality_summary(
        path_to_poscar: str | None = None,
        path_to_lobsterout: str | None = None,
        path_to_lobsterin: str | None = None,
        path_to_potcar: str | None = None,
        potcar_symbols: list | None = None,
        path_to_charge: str | None = None,
        path_to_bandoverlaps: str | None = None,
        path_to_doscar: str | None = None,
        path_to_vasprun: str | None = None,
        structure_obj: Structure | None = None,
        lobsterin_obj: Lobsterin | None = None,
        lobsterout_obj: Lobsterout | None = None,
        charge_obj: Charge | None = None,
        bandoverlaps_obj: Bandoverlaps | None = None,
        lobster_completedos_obj: LobsterCompleteDos | None = None,
        vasprun_obj: Vasprun | None = None,
        dos_comparison: bool = False,
        e_range: list = [-5, 0],
        n_bins: int | None = None,
        bva_comp: bool = False,
    ) -> dict:
        """
        Analyze LOBSTER calculation quality.

        :param path_to_poscar: path to structure file
        :param path_to_lobsterout: path to lobsterout file
        :param path_to_lobsterin: path to lobsterin file
        :param path_to_potcar: path to VASP potcar file
        :param potcar_symbols: list of potcar symbols from potcar file (can be used if no potcar available)
        :param path_to_charge: path to CHARGE.lobster file
        :param path_to_bandoverlaps: path to bandOverlaps.lobster file
        :param path_to_doscar: path to DOSCAR.lobster or DOSCAR.LSO.lobster file
        :param path_to_vasprun: path to vasprun.xml file
        :param structure_obj: pymatgen pymatgen.core.structure.Structure object
        :param lobsterin_obj: pymatgen.lobster.io.Lobsterin object
        :param lobsterout_obj: pymatgen lobster.io.Lobsterout object
        :param charge_obj: pymatgen lobster.io.Charge object
        :param bandoverlaps_obj: pymatgen lobster.io.BandOverlaps object
        :param lobster_completedos_obj: pymatgen.electronic_structure.dos.LobsterCompleteDos object
        :param vasprun_obj: pymatgen vasp.io.Vasprun object
        :param dos_comparison: will compare DOS from VASP and LOBSTER and return tanimoto index
        :param e_range: energy range for DOS comparisons
        :param n_bins: number of bins to discretize DOS for comparisons
        :param bva_comp: Compares LOBSTER charge signs with Bond valence charge signs

        Returns:
            A dict of summary of LOBSTER calculation quality by analyzing basis set used,
            charge spilling from lobsterout/ PDOS comparisons of VASP and LOBSTER /
            BVA charge comparisons

        """
        quality_dict = {}

        if path_to_potcar and not potcar_symbols and not path_to_vasprun and not vasprun_obj:
            potcar_names = Lobsterin._get_potcar_symbols(POTCAR_input=path_to_potcar)
        elif not path_to_potcar and not path_to_vasprun and not vasprun_obj and potcar_symbols:
            potcar_names = potcar_symbols
        elif path_to_vasprun and not vasprun_obj:
            vasprun = Vasprun(path_to_vasprun, parse_potcar_file=False, parse_eigen=False)
            potcar_names = [potcar.split(" ")[1] for potcar in vasprun.potcar_symbols]
        elif vasprun_obj and not path_to_vasprun:
            potcar_names = [potcar.split(" ")[1] for potcar in vasprun_obj.potcar_symbols]
        else:
            raise ValueError(
                "Please provide either path_to_potcar or list of "
                "potcar_symbols or path to vasprun.xml or vasprun object. "
                "Crucial to identify basis used for projections"
            )

        if path_to_poscar:
            struct = Structure.from_file(path_to_poscar)
        elif structure_obj:
            struct = structure_obj
        else:
            raise ValueError("Please provide path_to_poscar or structure_obj")

        ref_bases = Lobsterin.get_all_possible_basis_functions(structure=struct, potcar_symbols=potcar_names)

        if path_to_lobsterin:
            lobs_in = Lobsterin.from_file(path_to_lobsterin)
        elif lobsterin_obj:
            lobs_in = lobsterin_obj
        else:
            raise ValueError("Please provide path_to_lobsterin or lobsterin_obj")

        calc_basis = []
        for basis in lobs_in["basisfunctions"]:
            basis_sep = basis.split()[1:]
            basis_comb = " ".join(basis_sep)
            calc_basis.append(basis_comb)

        if calc_basis == list(ref_bases[0].values()):
            quality_dict["minimal_basis"] = True  # type: ignore
        else:
            quality_dict["minimal_basis"] = False  # type: ignore
            warnings.warn(
                "Consider rerunning the calc with the minimum basis as well. Choosing is "
                "larger basis set is recommended if you see a significant improvement of "
                "the charge spilling and material has non-zero band gap."
            )

        if path_to_lobsterout:
            lob_out = Lobsterout(path_to_lobsterout)
        elif lobsterout_obj:
            lob_out = lobsterout_obj
        else:
            raise ValueError("Please provide path_to_lobsterout or lobsterout_obj")

        quality_dict["charge_spilling"] = {
            "abs_charge_spilling": round((sum(lob_out.charge_spilling) / 2) * 100, 4),
            "abs_total_spilling": round((sum(lob_out.total_spilling) / 2) * 100, 4),
        }  # type: ignore

        if path_to_bandoverlaps is not None and not bandoverlaps_obj:
            band_overlaps = Bandoverlaps(filename=path_to_bandoverlaps) if Path(path_to_bandoverlaps).exists() else None
        elif path_to_bandoverlaps is None and bandoverlaps_obj:
            band_overlaps = bandoverlaps_obj
        else:
            band_overlaps = None

        if band_overlaps is not None:
            for line in lob_out.warning_lines:
                if "k-points could not be orthonormalized" in line:
                    total_kpoints = int(line.split(" ")[2])

            # store actual number of devations above pymatgen default limit of 0.1
            dev_val = []
            for dev in band_overlaps.max_deviation:
                if dev > 0.1:
                    dev_val.append(dev)

            quality_dict["band_overlaps_analysis"] = {  # type: ignore
                "file_exists": True,
                "limit_maxDeviation": 0.1,
                "has_good_quality_maxDeviation": band_overlaps.has_good_quality_maxDeviation(limit_maxDeviation=0.1),
                "max_deviation": round(max(band_overlaps.max_deviation), 4),
                "percent_kpoints_abv_limit": round((len(dev_val) / total_kpoints) * 100, 4),
            }

        else:
            quality_dict["band_overlaps_analysis"] = {  # type: ignore
                "file_exists": False,
                "limit_maxDeviation": None,
                "has_good_quality_maxDeviation": True,
                "max_deviation": None,
                "percent_kpoints_abv_limit": None,
            }

        if bva_comp:
            try:
                bond_valence = BVAnalyzer()

                bva_oxi = []
                if path_to_charge and not charge_obj:
                    lobs_charge = Charge(filename=path_to_charge)
                elif not path_to_charge and charge_obj:
                    lobs_charge = charge_obj
                else:
                    raise Exception("BVA comparison is requested, thus please provide path_to_charge or charge_obj")
                for i in bond_valence.get_valences(structure=struct):
                    if i >= 0:
                        bva_oxi.append("POS")
                    else:
                        bva_oxi.append("NEG")

                mull_oxi = []
                for i in lobs_charge.Mulliken:
                    if i >= 0:
                        mull_oxi.append("POS")
                    else:
                        mull_oxi.append("NEG")

                loew_oxi = []
                for i in lobs_charge.Loewdin:
                    if i >= 0:
                        loew_oxi.append("POS")
                    else:
                        loew_oxi.append("NEG")

                quality_dict["charge_comparisons"] = {}  # type: ignore
                if mull_oxi == bva_oxi:
                    quality_dict["charge_comparisons"]["bva_mulliken_agree"] = True  # type: ignore
                else:
                    quality_dict["charge_comparisons"]["bva_mulliken_agree"] = False  # type: ignore

                if mull_oxi == bva_oxi:
                    quality_dict["charge_comparisons"]["bva_loewdin_agree"] = True  # type: ignore
                else:
                    quality_dict["charge_comparisons"]["bva_loewdin_agree"] = False  # type: ignore

            except ValueError:
                quality_dict["charge_comparisons"] = {}  # type: ignore
                warnings.warn(
                    "Oxidation states from BVA analyzer cannot be determined. "
                    "Thus BVA charge comparison will be skipped"
                )
        if dos_comparison:
            if "LSO" not in str(path_to_doscar).split("."):
                warnings.warn(
                    "Consider using DOSCAR.LSO.lobster, as non LSO DOS from LOBSTER can have negative DOS values"
                )
            if path_to_doscar:
                doscar_lobster = Doscar(
                    doscar=path_to_doscar,
                    structure_file=path_to_poscar,
                    structure=structure_obj,
                )

                dos_lobster = doscar_lobster.completedos
            elif lobster_completedos_obj:
                dos_lobster = lobster_completedos_obj
            else:
                raise ValueError(
                    "Dos comparison is requested, so please provide either path_to_doscar or lobster_completedos_obj"
                )

            if path_to_vasprun:
                vasprun = Vasprun(path_to_vasprun, parse_potcar_file=False, parse_eigen=False)
            elif vasprun_obj:
                vasprun = vasprun_obj
            else:
                raise ValueError(
                    "Dos comparison is requested, so please provide either path to vasprun.xml or vasprun_obj"
                )
            dos_vasp = vasprun.complete_dos

            quality_dict["dos_comparisons"] = {}  # type: ignore

            for orb in dos_lobster.get_spd_dos():
                if e_range[0] >= min(dos_vasp.energies) and e_range[0] >= min(dos_lobster.energies):
                    min_e = e_range[0]
                else:
                    warnings.warn(
                        "Minimum energy range requested for DOS comparisons is not available "
                        "in VASP or LOBSTER calculation. Thus, setting min_e to -5 eV"
                    )
                    min_e = -5

                if e_range[-1] <= max(dos_vasp.energies) and e_range[-1] <= max(dos_lobster.energies):
                    max_e = e_range[-1]
                else:
                    warnings.warn(
                        "Maximum energy range requested for DOS comparisons is not available "
                        "in VASP or LOBSTER calculation. Thus, setting max_e to 0 eV"
                    )
                    max_e = 0

                if np.diff(dos_vasp.energies)[0] >= 0.1 or np.diff(dos_lobster.energies)[0] >= 0.1:
                    warnings.warn(
                        "Input DOS files have very few points in the energy interval and thus "
                        "comparisons will not be reliable. Please rerun the calculations with "
                        "higher number of DOS points. Set NEDOS and COHPSteps tags to >= 2000 in VASP and LOBSTER "
                        "calculations, respectively."
                    )

                if not n_bins:
                    n_bins = 56

                fp_lobster_orb = dos_lobster.get_dos_fp(
                    min_e=min_e,
                    max_e=max_e,
                    n_bins=n_bins,
                    normalize=True,
                    type=orb.name,
                )
                fp_vasp_orb = dos_vasp.get_dos_fp(
                    min_e=min_e,
                    max_e=max_e,
                    n_bins=n_bins,
                    normalize=True,
                    type=orb.name,
                )

                tani_orb = round(
                    dos_vasp.get_dos_fp_similarity(fp_lobster_orb, fp_vasp_orb, tanimoto=True),
                    4,
                )
                quality_dict["dos_comparisons"][f"tanimoto_orb_{orb.name}"] = tani_orb  # type: ignore

            fp_lobster = dos_lobster.get_dos_fp(
                min_e=min_e,
                max_e=max_e,
                n_bins=n_bins,
                normalize=True,
                type="summed_pdos",
            )
            fp_vasp = dos_vasp.get_dos_fp(
                min_e=min_e,
                max_e=max_e,
                n_bins=n_bins,
                normalize=True,
                type="summed_pdos",
            )

            tanimoto_summed = round(dos_vasp.get_dos_fp_similarity(fp_lobster, fp_vasp, tanimoto=True), 4)
            quality_dict["dos_comparisons"]["tanimoto_summed"] = tanimoto_summed  # type: ignore
            quality_dict["dos_comparisons"]["e_range"] = [min_e, max_e]  # type: ignore
            quality_dict["dos_comparisons"]["n_bins"] = n_bins  # type: ignore

        return quality_dict
