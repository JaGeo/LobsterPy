import collections

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Analysis:
    def __init__(self, path_to_poscar, path_to_icohplist, path_to_cohpcar, path_to_charge=None,
                 whichbonds="cation-anion", \
                 cutoff_icohp=0.1):

        self.path_to_poscar = path_to_poscar
        self.path_to_icohplist = path_to_icohplist
        self.path_to_cohpcar = path_to_cohpcar
        self.whichbonds = whichbonds
        self.cutoff_icohp = cutoff_icohp
        self.path_to_charge = path_to_charge
        self.setup_env()
        self.get_information_all_bonds()

        if path_to_charge is None:
            self.type_charge = "Valences"
        else:
            self.type_charge = "Mulliken"

        self.set_condensed_bonding_analysis()
        # TODO: get good data formats for maschine learning (graph with ICOHP values, for example?)!

    def setup_env(self):

        self.structure = Structure.from_file(self.path_to_poscar)
        sga = SpacegroupAnalyzer(structure=self.structure)
        symmetry_dataset = sga.get_symmetry_dataset()
        equivalent_sites = symmetry_dataset['equivalent_atoms']

        self.set_equivalent_sites = list(set(equivalent_sites))
        self.spg = symmetry_dataset['international']
        # What do I need for an automated analysis:
        #
        if self.whichbonds == "cation-anion":
            self.chemv = LobsterNeighbors(filename_ICOHP=self.path_to_icohplist,
                                          structure=Structure.from_file(self.path_to_poscar),
                                          additional_condition=1,
                                          perc_strength_ICOHP=self.cutoff_icohp,
                                          filename_CHARGE=self.path_to_charge,
                                          valences_from_charges=True
                                          )
        else:
            ValueError("only cation anion bonds implemented so far")

        # determine cations and anions
        self.lse = self.chemv.get_light_structure_environment(only_cation_environments=True)

    def get_information_all_bonds(self):

        self.set_inequivalent_cations = []
        self.set_coordinations_cations = []
        self.set_infos_cation_anion_bonds = []
        self.set_labels_cohps = []
        self.set_cohps = []
        # only_bonds_to
        # Write something to get all different kinds of anions?
        self.anion_types = self.chemv.get_anion_types()
        for ice, ce in enumerate(self.lse.coordination_environments):

            # only look at inequivalent sites (use of symmetry!)!
            if ice in self.set_equivalent_sites and ce[0]['ce_symbol'] is not None:
                self.set_inequivalent_cations.append(ice)
                ce = ce[0]['ce_symbol']
                self.set_coordinations_cations.append(ce)
                cation_anion_infos = self.chemv.get_info_icohps_to_neighbors([ice])
                self.set_infos_cation_anion_bonds.append(cation_anion_infos)

                aniontype_labels = []
                aniontype_cohps = []

                for anion in self.anion_types:
                    # print(anion)
                    labels, summedcohps = self.chemv.get_info_cohps_to_neighbors(self.path_to_cohpcar, [ice],
                                                                                 summed_spin_channels=True,
                                                                                 per_bond=False,
                                                                                 only_bonds_to=[str(anion)])
                    aniontype_labels.append(labels)
                    aniontype_cohps.append(summedcohps)

                self.set_labels_cohps.append(aniontype_labels)
                self.set_cohps.append(aniontype_cohps)

    @staticmethod
    def _get_strenghts_for_each_bond(pairs, strengths, cationname):
        dict_strenghts = {}
        for pair, strength in zip(pairs, strengths):
            # print(pair)
            # print(strength)
            new = [LobsterNeighbors._split_string(pair[0])[0], LobsterNeighbors._split_string(pair[1])[0]]

            new = Analysis._sort_name(new, cationname)
            # print(new2)
            string_here = new[0] + "-" + new[1]
            if not string_here in dict_strenghts:
                dict_strenghts[string_here] = []
            dict_strenghts[string_here].append(strength)
        return dict_strenghts

    @staticmethod
    def _sort_name(pair, cationname):

        new = []
        if pair[0] in cationname:
            new.append(pair[0])
            new.append(pair[1])

        elif pair[1] in cationname:
            new.append(pair[1])
            new.append(pair[0])
        return new

    def _get_antibdg_states(self, cohps, labels, namecation):
        dict_antibd = {}
        for label, cohp in zip(labels, cohps):
            new = label.split(' ')[2].split('-')
            sorted_new = self._sort_name(new, namecation)
            new_label = sorted_new[0] + '-' + sorted_new[1]
            antbd = cohp.has_antibnd_states_below_efermi(limit=0.1)
            if Spin.down in antbd:
                dict_antibd[new_label] = antbd[Spin.up] or antbd[Spin.down]
            else:
                dict_antibd[new_label] = antbd[Spin.up]

        return dict_antibd

    @staticmethod
    def _get_bond_dict(atoms, antbd):
        bond_dict = {}
        for key, item in atoms.items():
            bond_dict[key.split("-")[1]] = {"ICOHP_mean": str(round(np.mean(item), 2)),
                                            "has_antibdg_states_below_Efermi": antbd[key],
                                            "number_of_bonds": len(item)}

        return bond_dict

    @staticmethod
    def _get_bond_labels(atoms):
        all_labels = []
        for atomsnames in atoms:
            new = [LobsterNeighbors._split_string(atomsnames[0])[0], LobsterNeighbors._split_string(atomsnames[1])[0]]
            new.sort()
            # print(new2)
            string_here = new[0] + "-" + new[1]
            all_labels.append(string_here)
        count = collections.Counter(all_labels)
        return count

    def set_condensed_bonding_analysis(self):
        self.condensed_bonding_analysis = {}
        limit_icohps = self.chemv._get_limit_from_extremum(self.chemv.Icohpcollection, self.cutoff_icohp)

        formula = str(self.structure.composition.reduced_formula)
        number_considered_cations = len(self.set_inequivalent_cations)
        max_bond_lengths = max(self.chemv.Icohpcollection._list_length)
        charge_list = self.chemv.valences

        # information for each site: element, Mulliken charge, geometry, coordination environment, what kind of bonds,
        # sum of icohps for each of the bonds, bonding/antibonding interactions below e-fermi
        # nature of conduction/valence band?
        # information on Lobster run: what kind of basis, if possible ?

        site_dict = {}
        for ication, ce, cation_anion_infos, labels, cohps in zip(self.set_inequivalent_cations,
                                                                  self.set_coordinations_cations,
                                                                  self.set_infos_cation_anion_bonds,
                                                                  self.set_labels_cohps,
                                                                  self.set_cohps
                                                                  ):
            namecation = str(self.structure[ication].specie)

            # build something to get mean strengths for each type of bond
            mean_icohps = self._get_strenghts_for_each_bond(cation_anion_infos[4], cation_anion_infos[1], namecation)

            antbdg = self._get_antibdg_states(cohps, labels, namecation)
            bond_dict = self._get_bond_dict(mean_icohps, antbdg)

            site_dict[ication] = {"env": ce, "bonds": bond_dict, "cation": namecation, "charge": charge_list[ication]}

        self.condensed_bonding_analysis = {"formula": formula, "max_considered_bond_length": max_bond_lengths,
                                           "limit_icohp": limit_icohps, "number_of_considered_cations":
                                               number_considered_cations, "sites": site_dict,
                                           "type_charges": self.type_charge
                                           }

    def get_condensed_bonding_analysis(self):
        return self.condensed_bonding_analysis
