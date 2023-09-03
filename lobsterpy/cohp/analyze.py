# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines classes to analyze the COHPs automatically
"""
from __future__ import annotations

import warnings
from collections import Counter

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Analysis:
    """
    Analysis class of COHP data from Lobster

    Attributes:
        condensed_bonding_analysis: dict including a summary of the most important bonding properties
        final_dict_bonds: dict including information on ICOHPs per bond type
        final_dict_ions: dict including information on environments of cations
        chemenv: pymatgen.io.lobster.lobsterenv.LobsterNeighbors object
        lse: LightStructureEnvironment from pymatgen
        cutoff_icohp: Cutoff in percentage for evaluating neighbors based on ICOHP values.
         cutoff_icohp*max_icohp limits the number of considered environments
        anion_types: Set of Element objects from pymatgen
        list_equivalent_sites: list of site indices of sites that indicate which sites are equivalent
         e.g., [0 1 2 2 2] where site 0, 1, 2 indicate sites that are independent from each other
        path_to_charge: str that describes the path to CHARGE.lobster
        path_to_cohpcar: str that describes the path to COHPCAR.lobster
        path_to_icohplist: str that describes the path to ICOHPLIST.lobster
        path_to_poscar: str that describes path to POSCAR
        path_to_madelung: str that describes path to POSCAR
        seq_cohps: list of cohps
        seq_coord_ions: list of coodination environment strings for each cation
        set_equivalent_sites: set of inequivalent sites
        seq_ineq_ions: set of inequivalent cations/sites in the structure
        seq_infos_bonds (list): information on cation anion bonds (lists
            of pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo)
        spg: space group information
        structure: Structure object
        type_charge: which charges are considered here
        whichbonds: which bonds will be considered in analysis

    """

    def __init__(
        self,
        path_to_poscar: str,
        path_to_icohplist: str,
        path_to_cohpcar: str,
        path_to_charge: str | None = None,
        path_to_madelung: str | None = None,
        whichbonds: str = "cation-anion",
        cutoff_icohp: float = 0.1,
        summed_spins=True,
        type_charge=None,
        start=None,
    ):
        """
        This is a class to analyse bonding information automatically

        Args:
            path_to_poscar: path to POSCAR (e.g., "POSCAR")
            path_to_icohplist: path to ICOHPLIST.lobster (e.g., "ICOHPLIST.lobster")
            path_to_cohpcar: path to COHPCAR.lobster (e.g., "COHPCAR.lobster")
            path_to_charge: path to CHARGE.lobster (e.g., "CHARGE.lobster")
            path_to_madelung: path to MadelungEnergies.lobster (e.g., "MadelungEnergies.lobster")
            whichbonds: selects which kind of bonds are analyzed. "cation-anion" is the default
            cutoff_icohp: only bonds that are stronger than cutoff_icohp*strongest ICOHP will be considered
            summed_spins: if true, spins will be summed
            type_charge: If no path_to_charge is given, Valences will be used. Otherwise, Mulliken charges.
                        Löwdin charges cannot be selected at the moment.
            start: start energy for integration
        """
        self.start = start
        self.path_to_poscar = path_to_poscar
        self.path_to_icohplist = path_to_icohplist
        self.path_to_cohpcar = path_to_cohpcar
        self.whichbonds = whichbonds
        self.cutoff_icohp = cutoff_icohp
        self.path_to_charge = path_to_charge
        self.path_to_madelung = path_to_madelung
        self.setup_env()
        self.get_information_all_bonds(summed_spins=summed_spins)

        # This determines how cations and anions
        if path_to_charge is None:
            self.type_charge = "Valences"
        else:
            if type_charge is None:
                self.type_charge = "Mulliken"
            elif type_charge == "Mulliken":
                self.type_charge = "Mulliken"
            elif type_charge == "Löwdin":
                raise ValueError(
                    "Only Mulliken charges can be used here at the moment. Implementation will follow."
                )
            else:
                self.type_charge = "Valences"
                print(
                    "type_charge cannot be read! Please use Mulliken/Löwdin. Now, we will use valences"
                )

        self.set_condensed_bonding_analysis()
        self.set_summary_dicts()
        self.path_to_madelung = path_to_madelung

    def setup_env(self):
        """
        This method helps setting up the light structure environments based on COHPs

        Returns:
            None

        """
        self.structure = Structure.from_file(self.path_to_poscar)
        sga = SpacegroupAnalyzer(structure=self.structure)
        symmetry_dataset = sga.get_symmetry_dataset()
        equivalent_sites = symmetry_dataset["equivalent_atoms"]
        self.list_equivalent_sites = equivalent_sites
        self.set_equivalent_sites = list(set(equivalent_sites))
        self.spg = symmetry_dataset["international"]

        if self.whichbonds == "cation-anion":
            try:
                self.chemenv = LobsterNeighbors(
                    filename_ICOHP=self.path_to_icohplist,
                    structure=Structure.from_file(self.path_to_poscar),
                    additional_condition=1,
                    perc_strength_ICOHP=self.cutoff_icohp,
                    filename_CHARGE=self.path_to_charge,
                    valences_from_charges=True,
                    adapt_extremum_to_add_cond=True,
                )
            except ValueError as err:
                if (
                    str(err) == "min() arg is an empty sequence"
                    or str(err)
                    == "All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work"
                ):
                    raise ValueError(
                        "Consider switching to an analysis of all bonds and not only cation-anion bonds."
                        " It looks like no cations are detected."
                    )
                raise err
        elif self.whichbonds == "all":
            # raise ValueError("only cation anion bonds implemented so far")
            self.chemenv = LobsterNeighbors(
                filename_ICOHP=self.path_to_icohplist,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=0,
                perc_strength_ICOHP=self.cutoff_icohp,
                filename_CHARGE=self.path_to_charge,
                valences_from_charges=True,
                adapt_extremum_to_add_cond=True,
            )

        else:
            raise ValueError("only cation anion and all bonds implemented so far")

        # determine cations and anions
        try:
            if self.whichbonds == "cation-anion":
                self.lse = self.chemenv.get_light_structure_environment(
                    only_cation_environments=True
                )
            elif self.whichbonds == "all":
                self.lse = self.chemenv.get_light_structure_environment(
                    only_cation_environments=False
                )
        except ValueError:

            class Lse:
                """Test class when error was raised"""

                def __init__(self, chemenv, valences=None):
                    """
                    Test class when error was raised

                    Args:
                        chemenv (LobsterNeighbors): LobsterNeighbors object
                        valences: list of valences

                    """
                    if valences is None:
                        self.coordination_environments = [
                            [{"ce_symbol": str(len(coord))}] for coord in chemenv
                        ]
                    else:
                        self.coordination_environments = []

                        for val, coord in zip(valences, chemenv):
                            if val >= 0.0:
                                self.coordination_environments.append(
                                    [{"ce_symbol": str(len(coord))}]
                                )
                            else:
                                self.coordination_environments.append(
                                    [{"ce_symbol": None}]
                                )

            if self.whichbonds == "all":
                self.lse = Lse(self.chemenv.list_coords)
            elif self.whichbonds == "cation-anion":
                # make a new list
                self.lse = Lse(self.chemenv.list_coords, self.chemenv.valences)

    def get_information_all_bonds(self, summed_spins=True):
        """
        This method will gather all information on the bonds within the compound

        Returns:
            None

        """
        if self.whichbonds == "cation-anion":
            # this will only analyze cation anion bonds which simplifies the analysis
            self.seq_ineq_ions = []
            self.seq_coord_ions = []
            self.seq_infos_bonds = []
            self.seq_labels_cohps = []
            self.seq_cohps = []
            # only_bonds_to

            self.anion_types = self.chemenv.get_anion_types()
            for ice, ce in enumerate(self.lse.coordination_environments):
                # only look at inequivalent sites (use of symmetry to speed everything up!)!
                # only look at those cations that have cation-anion bonds
                if ice in self.set_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.seq_ineq_ions.append(ice)

                    self.seq_coord_ions.append(ce[0]["ce_symbol"])
                    self.seq_infos_bonds.append(
                        self.chemenv.get_info_icohps_to_neighbors([ice])
                    )

                    aniontype_labels = []
                    aniontype_cohps = []

                    # go through all anions in the structure!
                    for anion in self.anion_types:
                        # get labels and summed cohp objects
                        labels, summedcohps = self.chemenv.get_info_cohps_to_neighbors(
                            self.path_to_cohpcar,
                            [ice],
                            summed_spin_channels=summed_spins,
                            per_bond=False,
                            only_bonds_to=[str(anion)],
                        )

                        aniontype_labels.append(labels)
                        aniontype_cohps.append(summedcohps)

                    self.seq_labels_cohps.append(aniontype_labels)
                    self.seq_cohps.append(aniontype_cohps)

        elif self.whichbonds == "all":
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
                if ice in self.set_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.seq_ineq_ions.append(ice)
                    self.seq_coord_ions.append(ce[0]["ce_symbol"])
                    self.seq_infos_bonds.append(
                        self.chemenv.get_info_icohps_to_neighbors([ice])
                    )

                    type_labels = []
                    type_cohps = []

                    for element in self.elements:
                        # get labels and summed cohp objects
                        labels, summedcohps = self.chemenv.get_info_cohps_to_neighbors(
                            self.path_to_cohpcar,
                            [ice],
                            onlycation_isites=False,
                            summed_spin_channels=summed_spins,
                            per_bond=False,
                            only_bonds_to=[str(element)],
                        )

                        type_labels.append(labels)
                        type_cohps.append(summedcohps)

                    self.seq_labels_cohps.append(type_labels)
                    self.seq_cohps.append(type_cohps)

    @staticmethod
    def _get_strenghts_for_each_bond(pairs, strengths, nameion=None):
        """

        Args:
            pairs: list of list including labels for the atoms, e.g., [['O3', 'Cu1'], ['O3', 'Cu1']]
            strengths (list of float): list that gives the icohp strengths as a float, [-1.86287, -1.86288]
            nameion: string including the name of the cation in the list, e.g Cu1

        Returns:
            dict including inormation on icohps for each bond type, e.g.
            {'Yb-Sb': [-1.59769, -2.14723, -1.7925, -1.60773, -1.80149, -2.14335]}


        """
        dict_strenghts = {}

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
    def _sort_name(pair, nameion=None):
        """
        Will place the cation first in a list of name strings
        Args:
            pair: ["O","Cu"]
            nameion: "Cu"

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

    def _get_antibdg_states(self, cohps, labels, nameion=None, limit=0.01):
        """
        Will return a dictionary including information on antibonding states
        e.g., similar to: {'Cu-O': True, 'Cu-F': True}

        Args:
            cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
            labels: ['2 x Cu-O', '4 x Cu-F']
            nameion: string of the cation name, e.g. "Cu"
            limit: limit to detect antibonding states

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

    def _integrate_antbdstates_below_efermi_for_set_cohps(self, labels, cohps, nameion):
        """
        .. warning:: NEEDS MORE TESTS

        This method will return a dictionary including information on antibonding states
        important is however that only the energy range can be considered that has been computed
        (i.e., this might not be all)
        e.g., similar to: {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995},
        'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}

        Args:
            cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
            labels: ['2 x Cu-O', '4 x Cu-F']
            nameion: string of the cation name, e.g. "Cu"

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
                (
                    integral,
                    perc,
                    integral2,
                    perc2,
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

    def _integrate_antbdstates_below_efermi(self, cohp, start):
        """
        .. warning:: NEEDS MORE TESTS

        This integrates the whole COHP curve that has been computed.
        The energy range is very important.
        At present the energy range considered is dependent on COHPstartEnergy
        set during lobster runs. The bonding / antibonding intergral values are sensitive to this parameter.
        If COHPstartEnergy value does not cover entire range of VASP calculations then
        absolute value of ICOHP_sum might not be equivalent to (bonding- antibonding) integral values.

        Args:
            cohp: cohp object
            start: integration start energy in eV , eg start = -15

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

            This will integrate only bonding interactions of the COHP

            Args:
                y: COHP values
                x: Energy values

            Returns:
                integrated value of bonding interactions
            """
            y = np.asanyarray(y)
            x = np.asanyarray(x)

            bonding = trapezoid(y, x)

            return np.round(bonding, 2)

        def integrate_negative(y, x):
            """
            Will integrate only one side of the COHP
            Args:
                y: COHP values
                x: Energy values
            Returns:
                integrated value of antibonding interactions
            """
            y = np.asanyarray(y)
            x = np.asanyarray(x)
            antibonding = trapezoid(y, x)

            return np.round(antibonding, 2)

        # will integrate spin.up and spin.down only below efermi
        energies_corrected = cohp.energies - cohp.efermi
        if Spin.down in cohp.cohp:
            summedcohp = cohp.cohp[Spin.up] + cohp.cohp[Spin.down]
        else:
            summedcohp = cohp.cohp[Spin.up]

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

    @staticmethod
    def _get_bond_dict(
        bond_strength_dict, small_antbd_dict, nameion=None, large_antbd_dict=None
    ):
        """
        Will return a bond_dict including information for each site

        Args:
            bond_strength_dict (dict): dict with bond names as key and lists of bond strengths as items
            small_antbd_dict (dict): dict including if there are antibonding interactions, {'Yb-Sb': False}
            nameion (str): name of the cation, e.g. Yb
            large_antbdg_dict: will be implemented later


        Returns:
            dict including information on the anion (as label) and the ICOHPs in the item of the dict
            ICOHP_mean refers to the mean ICOHP in eV
            ICOHP_sum refers to the sum of the ICOHPs in eV
            has_antibdg_states_below_Efermi is True if there are antibonding interactions below Efermi
            "number_of_bonds" will count the numbers of bonds to the cation

            example:
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
                    "ICOHP_mean": str(round(np.mean(item), 2)),
                    "ICOHP_sum": str(round(np.sum(item), 2)),
                    "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                    "number_of_bonds": len(item),
                }
            else:
                bond_dict[key_here] = {
                    "ICOHP_mean": str(round(np.mean(item), 2)),
                    "ICOHP_sum": str(round(np.sum(item), 2)),
                    "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                    "number_of_bonds": len(item),
                    "perc_antibdg_states_below_Efermi": large_antbd_dict[key],
                }

        return bond_dict

    def set_condensed_bonding_analysis(self):
        """
        Sets a condensed version of the bonding analysis including a summary dictionary

        Returns:
            None

        """
        self.condensed_bonding_analysis = {}
        # which icohps are considered
        if self.whichbonds == "cation-anion":
            limit_icohps = self.chemenv._get_limit_from_extremum(
                self.chemenv.Icohpcollection,
                self.cutoff_icohp,
                adapt_extremum_to_add_cond=True,
                additional_condition=1,
            )
        elif self.whichbonds == "all":
            limit_icohps = self.chemenv._get_limit_from_extremum(
                self.chemenv.Icohpcollection,
                self.cutoff_icohp,
                adapt_extremum_to_add_cond=True,
                additional_condition=0,
            )
            # formula of the compound
        formula = str(self.structure.composition.reduced_formula)

        # how many inequivalent cations are in the structure
        if self.whichbonds == "cation-anion":
            number_considered_ions = len(self.seq_ineq_ions)
        elif self.whichbonds == "all":
            number_considered_ions = len(self.seq_ineq_ions)

        # what was the maximum bond lengths that was considered
        max_bond_lengths = max(self.chemenv.Icohpcollection._list_length)

        # what are the charges for the cations in the structure
        charge_list = self.chemenv.valences

        # dictionary including bonding information for each site
        site_dict = {}
        if self.whichbonds == "cation-anion":
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
                dict_antibonding = (
                    self._integrate_antbdstates_below_efermi_for_set_cohps(
                        labels, cohps, nameion=namecation
                    )
                )

                bond_dict = self._get_bond_dict(mean_icohps, antbdg, namecation)

                for k, v in bond_dict.items():
                    for k2, v2 in dict_antibonding.items():
                        if namecation == k2.split("-")[0] and k == k2.split("-")[1]:
                            v["bonding"] = v2["bonding"]
                            v["antibonding"] = v2["antibonding"]

                site_dict[ication] = {
                    "env": ce,
                    "bonds": bond_dict,
                    "ion": namecation,
                    "charge": charge_list[ication],
                    "relevant_bonds": cation_anion_infos[3],
                }
        elif self.whichbonds == "all":
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

                dict_antibonding = (
                    self._integrate_antbdstates_below_efermi_for_set_cohps(
                        labels, cohps, nameion
                    )
                )

                bond_dict = self._get_bond_dict(mean_icohps, antbdg, nameion=nameion)

                for k, v in bond_dict.items():
                    for k2, v2 in dict_antibonding.items():
                        if nameion == k2.split("-")[0] and k == k2.split("-")[1]:
                            v["bonding"] = v2["bonding"]
                            v["antibonding"] = v2["antibonding"]

                site_dict[iion] = {
                    "env": ce,
                    "bonds": bond_dict,
                    "ion": nameion,
                    "charge": charge_list[iion],
                    "relevant_bonds": bond_infos[3],
                }

        if self.path_to_madelung is None:
            if self.whichbonds == "cation-anion":
                # This sets the dictionary including the most important information on the compound
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    "limit_icohp": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                }
            elif self.whichbonds == "all":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    "limit_icohp": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                }
        else:
            from pymatgen.io.lobster import MadelungEnergies

            madelung = MadelungEnergies(self.path_to_madelung)
            if self.type_charge == "Mulliken":
                madelung_energy = madelung.madelungenergies_Mulliken
            elif self.type_charge == "Löwdin":
                madelung_energy = madelung.madelungenergies_Loewdin
            # This sets the dictionary including the most important information on the compound
            if self.whichbonds == "cation-anion":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    "limit_icohp": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                    "madelung_energy": madelung_energy,
                }
            elif self.whichbonds == "all":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    "limit_icohp": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                    "madelung_energy": madelung_energy,
                }

    def set_summary_dicts(self):
        """
        Sets summary dicts that can be used for correlations

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
        relevant_ion_ids = [
            isite for isite in self.list_equivalent_sites if isite in self.seq_ineq_ions
        ]

        # formula_units = self.structure.composition.num_atoms /
        # self.structure.composition.reduced_composition.num_atoms

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
                        "ICOHP_sum": float(properties["ICOHP_sum"]),
                        "has_antbdg": properties["has_antibdg_states_below_Efermi"],
                    }
                else:
                    final_dict_bonds[label]["number_of_bonds"] += int(
                        properties["number_of_bonds"]
                    )
                    final_dict_bonds[label]["ICOHP_sum"] += float(
                        properties["ICOHP_sum"]
                    )
                    final_dict_bonds[label]["has_antbdg"] = (
                        final_dict_bonds[label]["has_antbdg"]
                        or properties["has_antibdg_states_below_Efermi"]
                    )
        self.final_dict_bonds = {}
        for key, item in final_dict_bonds.items():
            self.final_dict_bonds[key] = {}
            self.final_dict_bonds[key]["ICOHP_mean"] = item["ICOHP_sum"] / (
                item["number_of_bonds"]
            )
            self.final_dict_bonds[key]["has_antbdg"] = item["has_antbdg"]

        # rework, add all environments!
        final_dict_ions = {}
        for key in relevant_ion_ids:
            if (
                self.condensed_bonding_analysis["sites"][key]["ion"]
                not in final_dict_ions
            ):
                final_dict_ions[
                    self.condensed_bonding_analysis["sites"][key]["ion"]
                ] = [self.condensed_bonding_analysis["sites"][key]["env"]]
            else:
                final_dict_ions[
                    self.condensed_bonding_analysis["sites"][key]["ion"]
                ].append(self.condensed_bonding_analysis["sites"][key]["env"])

        self.final_dict_ions = {}
        for key, item in final_dict_ions.items():
            self.final_dict_ions[key] = dict(Counter(item))
