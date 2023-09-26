# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines classes to analyze the COHPs automatically
"""
from __future__ import annotations
from pathlib import Path

import warnings
from collections import Counter
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.lobster import Lobsterin, Doscar, Lobsterout, Charge, Bandoverlaps
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
        are_cobis : bool indicating if file contains COBI/ICOBI data
        are_coops : bool indicating if file contains COOP/ICOOP data
        noise_cutoff : float that sets the lower limit of icohps or icoops or icobis considered
        set_cohps: list of cohps
        set_coordination_ions: list of coodination environment strings for each cation
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
        noise_cutoff: float = 0.1,
        summed_spins=True,
        are_cobis=False,
        are_coops=False,
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
            are_cobis : bool indicating if file contains COBI/ICOBI data
            are_coops : bool indicating if file contains COOP/ICOOP data
            noise_cutoff : float that sets the lower limit of icohps or icoops or icobis considered
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
        self.are_cobis = are_cobis
        self.are_coops = are_coops
        self.noise_cutoff = noise_cutoff
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
                    are_cobis=self.are_cobis,
                    are_coops=self.are_coops,
                    noise_cutoff=self.noise_cutoff,
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
                are_cobis=self.are_cobis,
                are_coops=self.are_coops,
                noise_cutoff=self.noise_cutoff,
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

    def _get_pop_type(self):
        """
        Convenience method to return the type of input population file

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
        bond_strength_dict,
        small_antbd_dict,
        nameion=None,
        large_antbd_dict=None,
        type_pop=None,
    ):
        """
        Will return a bond_dict including information for each site

        Args:
            bond_strength_dict (dict): dict with bond names as key and lists of bond strengths as items
            small_antbd_dict (dict): dict including if there are antibonding interactions, {'Yb-Sb': False}
            nameion (str): name of the cation, e.g. Yb
            large_antbdg_dict: will be implemented later
            type_pop: population type analyzed. eg. COHP


        Returns:
            Eg., if type_pop == 'COHP', will return
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
        # set population type
        type_pop = self._get_pop_type()
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

                bond_dict = self._get_bond_dict(
                    mean_icohps, antbdg, namecation, type_pop=type_pop
                )

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

                bond_dict = self._get_bond_dict(
                    mean_icohps, antbdg, nameion=nameion, type_pop=type_pop
                )

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
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                }
            elif self.whichbonds == "all":
                self.condensed_bonding_analysis = {
                    "formula": formula,
                    "max_considered_bond_length": max_bond_lengths,
                    f"limit_i{type_pop.lower()}": limit_icohps,
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
                    f"limit_i{type_pop.lower()}": limit_icohps,
                    "number_of_considered_ions": number_considered_ions,
                    "sites": site_dict,
                    "type_charges": self.type_charge,
                    "madelung_energy": madelung_energy,
                }
            elif self.whichbonds == "all":
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
                    final_dict_bonds[label]["number_of_bonds"] += int(
                        properties["number_of_bonds"]
                    )
                    final_dict_bonds[label][f"I{type_pop}_sum"] += float(
                        properties[f"I{type_pop}_sum"]
                    )
                    final_dict_bonds[label]["has_antbdg"] = (
                        final_dict_bonds[label]["has_antbdg"]
                        or properties["has_antibdg_states_below_Efermi"]
                    )
        self.final_dict_bonds = {}
        for key, item in final_dict_bonds.items():
            self.final_dict_bonds[key] = {}
            self.final_dict_bonds[key][f"I{type_pop}_mean"] = item[
                f"I{type_pop}_sum"
            ] / (item["number_of_bonds"])
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

    @staticmethod
    def get_lobster_calc_quality_summary(
        path_to_poscar: str,
        path_to_lobsterout: str,
        path_to_lobsterin: str,
        path_to_potcar: str | None = None,
        potcar_symbols: list | None = None,
        path_to_charge: str | None = None,
        path_to_bandoverlaps: str | None = None,
        path_to_doscar: str | None = None,
        path_to_vasprun: str | None = None,
        dos_comparison: bool = False,
        e_range: list = [-5, 0],
        n_bins: int | None = None,
        bva_comp: bool = False,
    ) -> dict:
        """
        This method will analyze LOBSTER calculation quality

        Args:
            path_to_poscar: path to structure file
            path_to_lobsterout: path to lobsterout file
            path_to_lobsterin: path to lobsterin file
            path_to_potcar: path to VASP potcar file
            potcar_symbols: list of potcar symbols from postcar file (can be used if no potcar available)
            path_to_charge: path to CHARGE.lobster file
            path_to_bandoverlaps: path to bandOverlaps.lobster file
            path_to_doscar: path to DOSCAR.lobster or DOSCAR.LSO.lobster file
            path_to_vasprun: path to vasprun.xml file
            dos_comparison: will compare DOS from VASP and LOBSTER and return tanimoto index
            e_range: energy range for DOS comparisons
            n_bins: number of bins to discretize DOS for comparisons
            bva_comp: Compares LOBSTER charge signs with Bond valence charge signs

        Returns:
            A dict of summary of LOBSTER calculation quality by analyzing basis set used,
            charge spilling from lobsterout/ PDOS comparisons of VASP and LOBSTER /
            BVA charge comparisons

        """
        quality_dict = {}

        if path_to_potcar and not potcar_symbols:
            potcar_names = Lobsterin._get_potcar_symbols(POTCAR_input=path_to_potcar)
        elif not path_to_potcar and potcar_symbols:
            potcar_names = potcar_symbols
        else:
            raise ValueError(
                "Please provide either path_to_potcar or list of "
                "potcar_symbols used for the calculations"
            )

        struct = Structure.from_file(path_to_poscar)

        ref_bases = Lobsterin.get_all_possible_basis_functions(
            structure=struct, potcar_symbols=potcar_names
        )

        lobs_in = Lobsterin.from_file(path_to_lobsterin)
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

        lob_out = Lobsterout(path_to_lobsterout)

        quality_dict["charge_spilling"] = {
            "abs_charge_spilling": round((sum(lob_out.charge_spilling) / 2) * 100, 4),
            "abs_total_spilling": round((sum(lob_out.total_spilling) / 2) * 100, 4),
        }  # type: ignore

        if path_to_bandoverlaps is not None:
            if Path(path_to_bandoverlaps).exists():  # type: ignore
                band_overlaps = Bandoverlaps(filename=path_to_bandoverlaps)
                for line in lob_out.warning_lines:
                    if "k-points could not be orthonormalized" in line:
                        total_kpoints = int(line.split(" ")[2])

                # store actual number of devations above pymatgen default limit of 0.1
                dev_val = []
                for dev in band_overlaps.max_deviation:
                    if dev > 0.1:
                        dev_val.append(dev)

                quality_dict["band_overlaps"] = {
                    "file_exists": True,
                    "limit_maxDeviation": 0.1,
                    "has_good_quality_maxDeviation": band_overlaps.has_good_quality_maxDeviation(
                        limit_maxDeviation=0.1
                    ),
                    "max_deviation": round(max(band_overlaps.max_deviation), 4),
                    "percent_kpoints_abv_limit": round(
                        (len(dev_val) / total_kpoints) * 100, 4
                    ),
                }  # type: ignore

            else:
                quality_dict["band_overlaps"] = {
                    "file_exists": False,
                    "limit_maxDeviation": None,
                    "has_good_quality_maxDeviation": True,
                    "max_deviation": None,
                    "percent_kpoints_abv_limit": None,
                }  # type: ignore

        if bva_comp:
            try:
                bond_valence = BVAnalyzer()

                bva_oxi = []
                lobs_charge = Charge(filename=path_to_charge)
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

                quality_dict["charges"] = {}  # type: ignore
                if mull_oxi == bva_oxi:
                    quality_dict["charges"]["BVA_Mulliken_agree"] = True  # type: ignore
                else:
                    quality_dict["charges"]["BVA_Mulliken_agree"] = False  # type: ignore

                if mull_oxi == bva_oxi:
                    quality_dict["charges"]["BVA_Loewdin_agree"] = True  # type: ignore
                else:
                    quality_dict["charges"]["BVA_Loewdin_agree"] = False  # type: ignore
            except ValueError:
                quality_dict["charges"] = {}  # type: ignore
                warnings.warn(
                    "Oxidation states from BVA analyzer cannot be determined. "
                    "Thus BVA charge comparison will be skipped"
                )
        if dos_comparison:
            if "LSO" not in str(path_to_doscar).split("."):
                warnings.warn(
                    "Consider using DOSCAR.LSO.lobster, as non LSO DOS from LOBSTER can have "
                    "negative DOS values"
                )
            doscar_lobster = Doscar(
                doscar=path_to_doscar,
                structure_file=path_to_poscar,
            )

            dos_lobster = doscar_lobster.completedos

            vasprun = Vasprun(path_to_vasprun)
            dos_vasp = vasprun.complete_dos

            quality_dict["dos_comparisons"] = {}  # type: ignore

            for orb in dos_lobster.get_spd_dos():
                if e_range[0] >= min(dos_vasp.energies) and e_range[0] >= min(
                    dos_lobster.energies
                ):
                    min_e = e_range[0]
                else:
                    warnings.warn(
                        "Mimimum energy range requested for DOS comparisons is not available "
                        "in VASP or LOBSTER calculation. Thus, setting min_e to -5 eV"
                    )
                    min_e = -5

                if e_range[-1] <= max(dos_vasp.energies) and e_range[-1] <= max(
                    dos_lobster.energies
                ):
                    max_e = e_range[-1]
                else:
                    warnings.warn(
                        "Maximum energy range requested for DOS comparisons is not available "
                        "in VASP or LOBSTER calculation. Thus, setting max_e to 0 eV"
                    )
                    max_e = 0

                if (
                    np.diff(dos_vasp.energies)[0] >= 0.1
                    or np.diff(dos_lobster.energies)[0] >= 0.1
                ):
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
                    dos_vasp.get_dos_fp_similarity(
                        fp_lobster_orb, fp_vasp_orb, tanimoto=True
                    ),
                    4,
                )
                quality_dict["dos_comparisons"][
                    "tanimoto_orb_{}".format(orb.name)
                ] = tani_orb  # type: ignore

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

            tanimoto_summed = round(
                dos_vasp.get_dos_fp_similarity(fp_lobster, fp_vasp, tanimoto=True), 4
            )
            quality_dict["dos_comparisons"]["tanimoto_summed"] = tanimoto_summed  # type: ignore
            quality_dict["dos_comparisons"]["e_range"] = [min_e, max_e]  # type: ignore
            quality_dict["dos_comparisons"]["n_bins"] = n_bins  # type: ignore

        return quality_dict
