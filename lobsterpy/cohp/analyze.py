# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines classes to analyze the COHPs automatically
"""

from collections import Counter
from typing import Optional

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Analysis:
    """
    Analysis class of COHP data from Lobster

    .. attribute:: condensed_bonding_analysis
        dict including a summary of the most important bonding properties

    .. attribute:: final_dict_bonds
        dict including information on ICOHPs per bond type

    .. attribute:: final_dict_ions
        dict including information on environments of cations

    .. attribute:: chemenv
        pymatgen.io.lobster.lobsterenv.LobsterNeighbors object

    .. attribute:: lse
        LightStructureEnvironment from pymatgen

    .. attribute:: cutoff_icohp
        Cutoff in percentage for evaluating neighbors based on ICOHP values.
        cutoff_icohp*max_icohp limits the number of considered environments

    .. attribute:: anion_types
        Set of Element objects from pymatgen

    .. attribute:: list_equivalent_sites
        list of site indices of sites that indicate which sites are equivalent
        e.g., [0 1 2 2 2] where site 0, 1, 2 indicate sites that are independent from each ther

    .. attribute:: path_to_charge
        str that describes the path to CHARGE.lobster

    .. attribute:: path_to_cohpcar
        str that describes the path to COHPCAR.lobster

    .. attribute:: path_to_icohplist
        str that describes the path to ICOHPLIST.lobster

    .. attribute:: path_to_poscar
        str that describes path to POSCAR

    .. attribute:: path_to_madelung
        str that describes path to POSCAR

    .. attribute:: set_cohps
        list of cohps

    .. attribute:: set_coordination_ions
        list of coodination environment strings for each cation

    .. attribute:: set_equivalent_sites
        set of inequivalent sites

    .. attribute:: set_inequivalent_ions
        set of inequivalent cations/sites in the structure

    .. attribute:: set_infos_bonds
        information on cation anion bonds

    .. attribute:: spg
        space group information

    .. attribute:: structure
        Structure object

    .. attribute:: type_charge
        which charges are considered here

    .. attribute:: whichbonds
        which bonds will be considered in analysis


    """

    def __init__(
        self,
        path_to_poscar: str,
        path_to_icohplist: str,
        path_to_cohpcar: str,
        path_to_charge: Optional[str] = None,
        path_to_madelung: Optional[str] = None,
        whichbonds: str = "cation-anion",
        cutoff_icohp: float = 0.1,
        summed_spins=True,
        type_charge=None,
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
                        Löwdin charges can be selected by using the keyword "Löwdin"
        """

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

        """
        self.structure = Structure.from_file(self.path_to_poscar)
        sga = SpacegroupAnalyzer(structure=self.structure)
        symmetry_dataset = sga.get_symmetry_dataset()
        equivalent_sites = symmetry_dataset["equivalent_atoms"]

        self.list_equivalent_sites = equivalent_sites
        self.set_equivalent_sites = list(set(equivalent_sites))
        self.spg = symmetry_dataset["international"]

        if self.whichbonds == "cation-anion":
            self.chemenv = LobsterNeighbors(
                filename_ICOHP=self.path_to_icohplist,
                structure=Structure.from_file(self.path_to_poscar),
                additional_condition=1,
                perc_strength_ICOHP=self.cutoff_icohp,
                filename_CHARGE=self.path_to_charge,
                valences_from_charges=True,
                adapt_extremum_to_add_cond=True,
            )

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
            raise ValueError("only cation anion bonds implemented so far")

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
                def __init__(self, chemenv):
                    self.coordination_environments = [
                        [{"ce_symbol": str(len(coord))}] for coord in chemenv
                    ]

            self.lse = Lse(self.chemenv.list_coords)

    def get_information_all_bonds(self, summed_spins=True):
        """
        This method will gather all information on the bonds within the compound
        Returns:

        """

        if self.whichbonds == "cation-anion":
            # this will only analyze cation anion bonds which simplifies the analysis
            self.set_inequivalent_ions = []
            self.set_coordination_ions = []
            self.set_infos_bonds = []
            self.set_labels_cohps = []
            self.set_cohps = []
            # only_bonds_to

            self.anion_types = self.chemenv.get_anion_types()
            for ice, ce in enumerate(self.lse.coordination_environments):
                # only look at inequivalent sites (use of symmetry to speed everything up!)!
                # only look at those cations that have cation-anion bonds
                if ice in self.set_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.set_inequivalent_ions.append(ice)
                    ce = ce[0]["ce_symbol"]
                    self.set_coordination_ions.append(ce)
                    cation_anion_infos = self.chemenv.get_info_icohps_to_neighbors(
                        [ice]
                    )
                    self.set_infos_bonds.append(cation_anion_infos)

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

                    self.set_labels_cohps.append(aniontype_labels)
                    self.set_cohps.append(aniontype_cohps)

        elif self.whichbonds == "all":
            # this will only analyze all bonds

            self.set_inequivalent_ions = []
            self.set_coordination_ions = []
            self.set_infos_bonds = []
            self.set_labels_cohps = []
            self.set_cohps = []
            # only_bonds_to
            self.elements = self.structure.composition.elements
            # self.anion_types = self.chemenv.get_anion_types()
            for ice, ce in enumerate(self.lse.coordination_environments):
                # only look at inequivalent sites (use of symmetry to speed everything up!)!
                # only look at those cations that have cation-anion bonds
                if ice in self.set_equivalent_sites and ce[0]["ce_symbol"] is not None:
                    self.set_inequivalent_ions.append(ice)
                    ce = ce[0]["ce_symbol"]
                    self.set_coordination_ions.append(ce)
                    bonds_infos = self.chemenv.get_info_icohps_to_neighbors([ice])
                    self.set_infos_bonds.append(bonds_infos)

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

                    self.set_labels_cohps.append(type_labels)
                    self.set_cohps.append(type_cohps)

    @staticmethod
    def _get_strenghts_for_each_bond(pairs, strengths, nameion=None):
        """

        Args:
            pairs: list of list including labels for the atoms, e.g., [['O3', 'Cu1'], ['O3', 'Cu1']]
            strengths (list of float): list that gives the icohp strenghts as a float, [-1.86287, -1.86288]
            nameion: string including the name of the cation in the list, e.g Cu1

        Returns: dict including inormation on icohps for each bond type, e.g.
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
                new = [
                    LobsterNeighbors._split_string(pair[0])[0],
                    LobsterNeighbors._split_string(pair[1])[0],
                ]
                new.sort()
                string_here = new[0] + "-" + new[1]

            if string_here not in dict_strenghts:
                dict_strenghts[string_here] = []
            dict_strenghts[string_here].append(strength)
        return dict_strenghts

    @staticmethod
    def _sort_name(pair, nameion=None):
        """
        will place the cation first in a list of name strings
        Args:
            pair: ["O","Cu"]
            nameion: "Cu"

        Returns: will return list of str, e.g. ["Cu", "O"]

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
        will return a dictionary including information on antibonding states
        e.g., similar to: {'Cu-O': True, 'Cu-F': True}
        Args:
            cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
            labels: ['2 x Cu-O', '4 x Cu-F']
            nameion: string of the cation name, e.g. "Cu"
            limit: limit to detect antibonding states

        Returns:    dict including in formation on whether antibonding interactions exist,
                    e.g., {'Cu-O': True, 'Cu-F': True}


        """

        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            # print(labels)
            if label is not None:
                if nameion is not None:
                    new = label.split(" ")[2].split("-")
                    sorted_new = self._sort_name(new, nameion)
                    new_label = sorted_new[0] + "-" + sorted_new[1]
                else:
                    new = label.split(" ")[2].split("-")
                    sorted_new = new.copy()
                    sorted_new.sort()
                    new_label = sorted_new[0] + "-" + sorted_new[1]

                antbd = cohp.has_antibnd_states_below_efermi(limit=limit)
                if Spin.down in antbd:
                    dict_antibd[new_label] = antbd[Spin.up] or antbd[Spin.down]
                else:
                    dict_antibd[new_label] = antbd[Spin.up]

        return dict_antibd

    def _integrate_antbdstates_below_efermi_for_set_cohps(self, labels, cohps, nameion):
        """
        WARNING: NEEDS MORE TESTS
        will return a dictionary including information on antibonding states
        important is however that only the energy range can be considered that has been computed
        (i.e., this might not be all)
        e.g., similar to: {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995}, 'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}
        Args:
            cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
            labels: ['2 x Cu-O', '4 x Cu-F']
            nameion: string of the cation name, e.g. "Cu"

        Returns:    dict including in formation on whether antibonding interactions exist,
                    e.g., {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995}, 'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}}
        """

        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            if label is not None:
                new = label.split(" ")[2].split("-")
                sorted_new = self._sort_name(new, nameion)
                new_label = sorted_new[0] + "-" + sorted_new[1]
                integral, perc = self._integrate_antbdstates_below_efermi(cohp, -2)
                dict_antibd[new_label] = {"integral": integral, "perc": perc}

        return dict_antibd

    def _integrate_antbdstates_below_efermi(self, cohp, start=-30):

        """
        WARNING: NEEDS MORE TESTS
        This integrates the whole COHP curve that has been computed. The energy range might be very important
        Args:
            cohp: cohp object
            start: where does the integration start

        Returns:
            absolute value of antibonding interactions, percentage value of antibonding interaction
        """

        # This integrates the whole COHP curve that has been computed. Just be aware that you might be
        # neglecting some low-lying interactions due to the energy range

        def abstrapz_positive(y, x=None, dx=0.001):
            """
            will integrate only one side of the COHP
            Args:
                y: Energy values
                x: COHP values
                dx: how fine should the integration steps be
            Returns:
                integrated value
            """

            y = np.asanyarray(y)
            if x is None:
                d = dx
            else:
                x = np.asanyarray(x)
                d = np.diff(x)
            ret = d * (y[1:] + y[:-1]) / 2.0
            return ret[ret > 0.0].sum()  # The important line

        def abstrapz_negative(y, x=None, dx=0.001):
            """
            will integrate only one side of the COHP
            Args:
                y: Energy values
                x: COHP values
                dx: how fine should the integration steps be

            Returns:
                integrated value
            """
            y = np.asanyarray(y)
            if x is None:
                d = dx
            else:
                x = np.asanyarray(x)
                d = np.diff(x)
            ret = d * (y[1:] + y[:-1]) / 2.0
            return ret[ret < 0.0].sum()  # The important line

        from scipy.interpolate import InterpolatedUnivariateSpline

        # will integrate spin.up and spin.down only below efermi and only below curve
        energies_corrected = cohp.energies - cohp.efermi
        if Spin.down in cohp.cohp:
            summedcohp = cohp.cohp[Spin.up] + cohp.cohp[Spin.down]
        else:
            summedcohp = cohp.cohp[Spin.up]

        spl = InterpolatedUnivariateSpline(energies_corrected, summedcohp, ext=0)
        # integrate only below curve
        # TODO: has to be tested again!
        integrate = abstrapz_negative(
            [spl(energy) for energy in energies_corrected if start <= energy <= 0],
            [energy for energy in energies_corrected if start <= energy <= 0],
        )

        integrate2 = abstrapz_positive(
            [spl(energy) for energy in energies_corrected if start <= energy <= 0],
            [energy for energy in energies_corrected if start <= energy <= 0],
        )

        return integrate2, abs(integrate2) / (abs(integrate2) + abs(integrate))

    @staticmethod
    def _get_bond_dict(
        bond_strength_dict, small_antbd_dict, nameion=None, large_antbd_dict=None
    ):
        """
        will return a bond_dict incluing information for each site
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
        sets a condensed version of the bonding analysis including a summary dictionary
        Returns:

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
            number_considered_ions = len(self.set_inequivalent_ions)
        elif self.whichbonds == "all":
            number_considered_ions = len(self.set_inequivalent_ions)

        # what was the maximum bond lengths that was considered
        max_bond_lengths = max(self.chemenv.Icohpcollection._list_length)

        # what are the charges for the cations in the structure
        charge_list = self.chemenv.valences

        # dictionary including bonding information for each site
        site_dict = {}
        if self.whichbonds == "cation-anion":
            for ication, ce, cation_anion_infos, labels, cohps in zip(
                self.set_inequivalent_ions,
                self.set_coordination_ions,
                self.set_infos_bonds,
                self.set_labels_cohps,
                self.set_cohps,
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
                # dict_antibonding = self._integrate_antbdstates_below_efermi_for_set_cohps(labels, cohps, nameion)
                # dict_antibonding,
                bond_dict = self._get_bond_dict(mean_icohps, antbdg, namecation)

                site_dict[ication] = {
                    "env": ce,
                    "bonds": bond_dict,
                    "ion": namecation,
                    "charge": charge_list[ication],
                    "relevant_bonds": cation_anion_infos[3],
                }
        elif self.whichbonds == "all":
            for iion, ce, bond_infos, labels, cohps in zip(
                self.set_inequivalent_ions,
                self.set_coordination_ions,
                self.set_infos_bonds,
                self.set_labels_cohps,
                self.set_cohps,
            ):
                nameion = str(self.structure[iion].specie)

                # This will compute the mean strengths of ICOHPs
                mean_icohps = self._get_strenghts_for_each_bond(
                    pairs=bond_infos[4], strengths=bond_infos[1], nameion=None
                )
                # pairs, strengths, nameion
                # will collect if there are antibonding states present
                antbdg = self._get_antibdg_states(cohps, labels, nameion=None)

                # dict_antibonding = self._integrate_antbdstates_below_efermi_for_set_cohps(labels, cohps, nameion)
                # dict_antibonding,
                bond_dict = self._get_bond_dict(mean_icohps, antbdg, nameion=nameion)

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
        sets summary dicts that can be used for correlations

        bond_dict that includes information on each bond

        "has_antbd" tells if there are antbonding states
        "ICOHP_mean" shows the mean of all ICOHPs in EV

        {'Yb-Sb': { 'has_antbdg': False, 'ICOHP_mean': -1.7448},
        'Mn-Sb': { 'has_antbdg': True, 'ICOHP_mean': -1.525}}

        a cation dict that includes all different coordination environments and counts for them
        {'Na': {'T:4': 4, 'A:2': 4}, 'Si': {'T:6': 4, 'PP:6': 4}}
        Returns:




        """
        relevant_ion_ids = [
            isite
            for isite in self.list_equivalent_sites
            if isite in self.set_inequivalent_ions
        ]

        # formula_units = self.structure.composition.num_atoms / self.structure.composition.reduced_composition.num_atoms

        final_dict_bonds = {}
        for key in relevant_ion_ids:

            item = self.condensed_bonding_analysis["sites"][key]
            for type, properties in item["bonds"].items():
                label_list = [item["ion"], str(type)]
                new_label = label_list.copy()
                new_label.sort()
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
