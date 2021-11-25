from typing import Optional

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# TODO: reduce the number of attributes

class Analysis:
    """
    Analysis class of COHP data from Lobster

    .. attribute:: condensed_bonding_analysis
        dict including a summary of the most important bonding properties

    .. attribute:: final_dict_bonds
        dict including information on ICOHPs per bond type

    .. attribute:: final_dict_cations
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

    .. attribute:: set_cohps
        list of cohps

    .. attribute:: set_coordinations_cations
        list of coodination environment strings for each cation

    .. attribute:: set_equivalent_sites
        set of inequivalent sites

    .. attribute:: set_inequivalent_cations
        set of inequivalent cations in the structure


    .. attribute:: set_infos_cation_anion_bonds
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

    def __init__(self, path_to_poscar: str, path_to_icohplist: str, path_to_cohpcar: str,
                 path_to_charge: Optional[str] = None,
                 whichbonds: str = "cation-anion", \
                 cutoff_icohp: float = 0.1):
        """
        This is a class to analyse bonding information automatically
        Args:
            path_to_poscar: path to POSCAR (e.g., "POSCAR")
            path_to_icohplist: path to ICOHPLIST.lobster (e.g., "ICOHPLIST.lobster")
            path_to_cohpcar: path to COHPCAR.lobster (e.g., "COHPCAR.lobster")
            path_to_charge: path to CHARGE.lobster (e.g., "CHARGE.lobster")
            whichbonds: selects which kind of bonds are analyzed. "cation-anion" is the default
            cutoff_icohp: only bonds that are stronger than cutoff_icohp*strongest ICOHP will be considered
        """

        self.path_to_poscar = path_to_poscar
        self.path_to_icohplist = path_to_icohplist
        self.path_to_cohpcar = path_to_cohpcar
        self.whichbonds = whichbonds
        self.cutoff_icohp = cutoff_icohp
        self.path_to_charge = path_to_charge
        self.setup_env()
        self.get_information_all_bonds()

        # This determines how cations and anions
        if path_to_charge is None:
            self.type_charge = "Valences"
        else:
            self.type_charge = "Mulliken"

        self.set_condensed_bonding_analysis()
        self.set_summary_dicts()

    def setup_env(self):
        """
        This method helps setting up the light structure environments based on COHPs
        Returns:

        """
        self.structure = Structure.from_file(self.path_to_poscar)
        sga = SpacegroupAnalyzer(structure=self.structure)
        symmetry_dataset = sga.get_symmetry_dataset()
        equivalent_sites = symmetry_dataset['equivalent_atoms']

        self.list_equivalent_sites = equivalent_sites
        self.set_equivalent_sites = list(set(equivalent_sites))
        self.spg = symmetry_dataset['international']
        # What do I need for an automated analysis:
        #
        if self.whichbonds == "cation-anion":
            self.chemenv = LobsterNeighbors(filename_ICOHP=self.path_to_icohplist,
                                            structure=Structure.from_file(self.path_to_poscar),
                                            additional_condition=1,
                                            perc_strength_ICOHP=self.cutoff_icohp,
                                            filename_CHARGE=self.path_to_charge,
                                            valences_from_charges=True,
                                            adapt_extremum_to_add_cond=True,
                                            )


        else:
            raise ValueError("only cation anion bonds implemented so far")

        # determine cations and anions
        self.lse = self.chemenv.get_light_structure_environment(only_cation_environments=True)

    def get_information_all_bonds(self):
        """
        This method will gather all information on the bonds within the compound
        Returns:

        """
        self.set_inequivalent_cations = []
        self.set_coordinations_cations = []
        self.set_infos_cation_anion_bonds = []
        self.set_labels_cohps = []
        self.set_cohps = []
        # only_bonds_to

        self.anion_types = self.chemenv.get_anion_types()
        for ice, ce in enumerate(self.lse.coordination_environments):
            # only look at inequivalent sites (use of symmetry to speed everything up!)!
            # only look at those cations that have cation-anion bonds
            if ice in self.set_equivalent_sites and ce[0]['ce_symbol'] is not None:
                self.set_inequivalent_cations.append(ice)
                ce = ce[0]['ce_symbol']
                self.set_coordinations_cations.append(ce)
                cation_anion_infos = self.chemenv.get_info_icohps_to_neighbors([ice])
                self.set_infos_cation_anion_bonds.append(cation_anion_infos)

                aniontype_labels = []
                aniontype_cohps = []

                # go through all anions in the structure
                # TODO: add options for other bonds as well
                if self.whichbonds == "cation-anion":
                    for anion in self.anion_types:
                        # get labels and summed cohp objects
                        labels, summedcohps = self.chemenv.get_info_cohps_to_neighbors(self.path_to_cohpcar, [ice],
                                                                                       summed_spin_channels=True,
                                                                                       per_bond=False,
                                                                                       only_bonds_to=[str(anion)])

                        aniontype_labels.append(labels)
                        aniontype_cohps.append(summedcohps)

                self.set_labels_cohps.append(aniontype_labels)
                self.set_cohps.append(aniontype_cohps)

    @staticmethod
    def _get_strenghts_for_each_bond(pairs, strengths, cationname):
        """

        Args:
            pairs: list of list including labels for the atoms, e.g., [['O3', 'Cu1'], ['O3', 'Cu1']]
            strengths (list of float): list that gives the icohp strenghts as a float, [-1.86287, -1.86288]
            cationname: string including the name of the cation in the list, e.g Cu1

        Returns: dict including inormation on icohps for each bond type, e.g.
        {'Yb-Sb': [-1.59769, -2.14723, -1.7925, -1.60773, -1.80149, -2.14335]}


        """
        dict_strenghts = {}

        for pair, strength in zip(pairs, strengths):
            new = [LobsterNeighbors._split_string(pair[0])[0], LobsterNeighbors._split_string(pair[1])[0]]
            new = Analysis._sort_name(new, cationname)
            string_here = new[0] + "-" + new[1]
            if not string_here in dict_strenghts:
                dict_strenghts[string_here] = []
            dict_strenghts[string_here].append(strength)
        return dict_strenghts

    @staticmethod
    def _sort_name(pair, cationname):
        """
        will place the cation first in a list of name strings
        Args:
            pair: ["O","Cu"]
            cationname: "Cu"

        Returns: will return list of str, e.g. ["Cu", "O"]

        """

        new = []
        if pair[0] == cationname:
            new.append(pair[0])
            new.append(pair[1])

        elif pair[1] == cationname:
            new.append(pair[1])
            new.append(pair[0])
        return new

    def _get_antibdg_states(self, cohps, labels, namecation, limit=0.01):
        """
        will return a dictionary including information on antibonding states
        e.g., similar to: {'Cu-O': True, 'Cu-F': True}
        Args:
            cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
            labels: ['2 x Cu-O', '4 x Cu-F']
            namecation: string of the cation name, e.g. "Cu"
            limit: limit to detect antibonding states

        Returns:    dict including in formation on whether antibonding interactions exist,
                    e.g., {'Cu-O': True, 'Cu-F': True}


        """

        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            if label is not None:
                new = label.split(' ')[2].split('-')
                sorted_new = self._sort_name(new, namecation)
                new_label = sorted_new[0] + '-' + sorted_new[1]
                antbd = cohp.has_antibnd_states_below_efermi(limit=limit)
                if Spin.down in antbd:
                    dict_antibd[new_label] = antbd[Spin.up] or antbd[Spin.down]
                else:
                    dict_antibd[new_label] = antbd[Spin.up]

        return dict_antibd

    def _integrate_antbdstates_below_efermi_for_set_cohps(self, labels, cohps, namecation):
        """
                will return a dictionary including information on antibonding states
                important is however that only the energy range can be considered that has been computed
                (i.e., this might not be all)
                e.g., similar to: {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995}, 'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}
                Args:
                    cohps: list of pymatgen.electronic_structure.cohp.Cohp ojbects
                    labels: ['2 x Cu-O', '4 x Cu-F']
                    namecation: string of the cation name, e.g. "Cu"

                Returns:    dict including in formation on whether antibonding interactions exist,
                            e.g., {'Cu-O': {'integral': 4.24374775705, 'perc': 5.7437713186999995}, 'Cu-F': {'integral': 3.07098300965, 'perc': 4.25800841445}}}
        """

        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            if label is not None:
                new = label.split(' ')[2].split('-')
                sorted_new = self._sort_name(new, namecation)
                new_label = sorted_new[0] + '-' + sorted_new[1]
                integral, perc = self._integrate_antbdstates_below_efermi(cohp, -2)
                dict_antibd[new_label] = {"integral": integral, "perc": perc}

        return dict_antibd

    def _integrate_antbdstates_below_efermi(self, cohp, start=-30):

        """
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
            ret = (d * (y[1:] + y[:-1]) / 2.0)
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
            ret = (d * (y[1:] + y[:-1]) / 2.0)
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
        integrate = abstrapz_negative([spl(energy) for energy in energies_corrected if start <= energy <= 0],
                                      [energy for energy in energies_corrected if start <= energy <= 0])

        integrate2 = abstrapz_positive([spl(energy) for energy in energies_corrected if start <= energy <= 0],
                                       [energy for energy in energies_corrected if start <= energy <= 0])

        return integrate2, abs(integrate2) / (abs(integrate2) + abs(integrate))

    @staticmethod
    def _get_bond_dict(bond_strength_dict, small_antbd_dict, namecation, large_antbd_dict=None):
        """
        will return a bond_dict incluing information for each site
        Args:
            bond_strength_dict (dict): dict with bond names as key and lists of bond strengths as items
            small_antbd_dict (dict): dict including if there are antibonding interactions, {'Yb-Sb': False}
            namecation (str): name of the cation, e.g. Yb
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

            a = key.split("-")[0]
            b = key.split("-")[1]
            if a == namecation:
                key_here = b
            elif b == namecation:
                key_here = a

            if large_antbd_dict is None:
                bond_dict[key_here] = {"ICOHP_mean": str(round(np.mean(item), 2)),
                                       "ICOHP_sum": str(round(np.sum(item), 2)),
                                       "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                                       "number_of_bonds": len(item),
                                       }
            else:
                bond_dict[key_here] = {"ICOHP_mean": str(round(np.mean(item), 2)),
                                       "ICOHP_sum": str(round(np.sum(item), 2)),
                                       "has_antibdg_states_below_Efermi": small_antbd_dict[key],
                                       "number_of_bonds": len(item),
                                       "perc_antibdg_states_below_Efermi": large_antbd_dict[key]}

        return bond_dict

    def set_condensed_bonding_analysis(self):
        """
        sets a condensed version of the bonding analysis including a summary dictionary
        Returns:

        """

        self.condensed_bonding_analysis = {}
        # which icohps are considered
        if self.whichbonds == "cation-anion":
            limit_icohps = self.chemenv._get_limit_from_extremum(self.chemenv.Icohpcollection, self.cutoff_icohp, adapt_extremum_to_add_cond=True, additional_condition=1)

        # formula of the compound
        formula = str(self.structure.composition.reduced_formula)

        # how many inequivalent cations are in the structure
        number_considered_cations = len(self.set_inequivalent_cations)

        # what was the maximum bond lengths that was considered
        max_bond_lengths = max(self.chemenv.Icohpcollection._list_length)

        # what are the charges for the cations in the structure
        charge_list = self.chemenv.valences

        # dictionary including bonding information for each site
        site_dict = {}
        for ication, ce, cation_anion_infos, labels, cohps in zip(self.set_inequivalent_cations,
                                                                  self.set_coordinations_cations,
                                                                  self.set_infos_cation_anion_bonds,
                                                                  self.set_labels_cohps,
                                                                  self.set_cohps
                                                                  ):
            namecation = str(self.structure[ication].specie)

            # This will compute the mean strengths of ICOHPs
            mean_icohps = self._get_strenghts_for_each_bond(pairs=cation_anion_infos[4],
                                                            strengths=cation_anion_infos[1], cationname=namecation)
            # pairs, strengths, cationname
            # will collect if there are antibonding states present
            antbdg = self._get_antibdg_states(cohps, labels, namecation)

            # dict_antibonding = self._integrate_antbdstates_below_efermi_for_set_cohps(labels, cohps, namecation)
            # dict_antibonding,
            bond_dict = self._get_bond_dict(mean_icohps, antbdg, namecation)

            site_dict[ication] = {"env": ce, "bonds": bond_dict, "cation": namecation, "charge": charge_list[ication],
                                  "relevant_bonds": cation_anion_infos[3]}

        # This sets the dictionary including the most important information on the compound
        self.condensed_bonding_analysis = {"formula": formula, "max_considered_bond_length": max_bond_lengths,
                                           "limit_icohp": limit_icohps, "number_of_considered_cations":
                                               number_considered_cations, "sites": site_dict,
                                           "type_charges": self.type_charge
                                           }

    def set_summary_dicts(self):
        """
        sets summary dicts that can be used for correlations

        bond_dict that includes information on each bond
        "number_of_bonds" counts the number of bonds
        "ICOHP_sum" shows the sum of all ICOHPs in eV
        "has_antbd" tells if there are antbonding states
        "ICOHP_mean" shows the mean of all ICOHPs in EV

        {'Yb-Sb': {'number_of_bonds': 25, 'ICOHP_sum': -43.62, 'has_antbdg': False, 'ICOHP_mean': -1.7448},
        'Mn-Sb': {'number_of_bonds': 4, 'ICOHP_sum': -6.1, 'has_antbdg': True, 'ICOHP_mean': -1.525}}

        a cation dict that includes all different coordination environments for symmetry inequivalent cations
        {'Yb': ['O:6', 'O:6', 'O:6', 'PB:7'], 'Mn': ['T:4']}
        Returns:




        """

        relevant_cation_ids = set([isite for isite in
                                   self.list_equivalent_sites if
                                   isite in self.set_inequivalent_cations])

        final_dict_bonds = {}
        for key in relevant_cation_ids:
            item = self.condensed_bonding_analysis["sites"][key]
            for type, properties in item["bonds"].items():
                label = item["cation"] + '-' + str(type)
                if not label in final_dict_bonds:
                    final_dict_bonds[label] = {"number_of_bonds": int(properties["number_of_bonds"]),
                                               "ICOHP_sum": float(properties["ICOHP_sum"]),
                                               "has_antbdg": properties["has_antibdg_states_below_Efermi"]}
                else:
                    final_dict_bonds[label]["number_of_bonds"] += int(properties["number_of_bonds"])
                    final_dict_bonds[label]["ICOHP_sum"] += float(properties["ICOHP_sum"])
                    final_dict_bonds[label]["has_antbdg"] = final_dict_bonds[label]["has_antbdg"] or properties[
                        "has_antibdg_states_below_Efermi"]

        for key, item in final_dict_bonds.items():
            final_dict_bonds[key]["ICOHP_mean"] = item["ICOHP_sum"] / item["number_of_bonds"]
            final_dict_bonds[key]["number_of_bonds"] = item["number_of_bonds"]

        self.final_dict_bonds = final_dict_bonds

        final_dict_cations = {}
        for key, item in self.condensed_bonding_analysis["sites"].items():

            if item["cation"] not in final_dict_cations:
                final_dict_cations[item["cation"]] = [item["env"]]
            else:
                final_dict_cations[item["cation"]].append(item["env"])

        self.final_dict_cations = final_dict_cations
