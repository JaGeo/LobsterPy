# TODO: write classes that allow to decribe the ICOHP situation in a crystal structure
# TODO: clean up!

# TODO: add a sentence which ICOHPs will be still considered

# TODO: einen Teil auslagern, um dann nur noch die Beschreibung hier zu haben?


import collections

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Description:
    """
    will write a description text of all relevant bonds. It will analyse all coordination environments.
    Furthermore, it will analyse which kind of bonds contribute to the valence band close to the Fermi level

    """

    def __init__(self, path_to_poscar, path_to_icohplist, path_to_cohpcar):

        # TODO: einen Teil auslagern in "Analysis" und dann als Startpunkt für "Description/Plot" benutzen
        # wichtig wäre auch, solche Datenformate zu bauen, die dann für maschinelles Lernen verwendet werden können.

        chemv, lse, set_equivalent_sites, structure = self.setup_env(path_to_icohplist, path_to_poscar)

        set_cohps, set_coordinations_cations, set_inequivalent_cations, set_infos_cation_anion_bonds, \
        set_labels_cohps, anion_types = self.get_information_all_bonds(
            chemv, lse, path_to_cohpcar, set_equivalent_sites)

        # TODO: wie nummeriere ich die Si durch
        self.set_description(set_cohps, set_coordinations_cations, set_inequivalent_cations,
                             set_infos_cation_anion_bonds, set_labels_cohps, structure)

    def plot_cohps(self, set_cohps, set_inequivalent_cations,
                   set_labels_cohps, structure):

        # TODO: check if self. can be used to simplify this plotting here
        for ication, labels, cohps in zip(set_inequivalent_cations, set_labels_cohps,
                                          set_cohps):

            namecation = str(structure[ication].specie)

            cp = CohpPlotter()
            for label, cohp in zip(labels, cohps):
                cp.add_cohp(namecation + str(ication) + ': ' + label, cohp)
            plot = cp.get_plot(integrated=False)
            plot.ylim([-4, 2])
            plot.show()

    @staticmethod
    def _get_strenghts_for_each_bond(pairs, strenghts, cationname):
        dict_strenghts = {}
        for pair, strength in zip(pairs, strenghts):
            # print(pair)
            # print(strength)
            new = [LobsterNeighbors._split_string(pair[0])[0], LobsterNeighbors._split_string(pair[1])[0]]

            new = Description._sort_name(new, cationname)
            # print(new2)
            string_here = new[0] + "-" + new[1]
            if not string_here in dict_strenghts:
                dict_strenghts[string_here] = []
            dict_strenghts[string_here].append(strength)
        return dict_strenghts

    @staticmethod
    def _sort_name(pair, cationname):

        new = ["", ""]

        if pair[0] in cationname:
            new[0] = pair[0]
            new[1] = pair[1]
        elif pair[1] in cationname:
            new[1] = pair[0]
            new[0] = pair[1]
        else:
            pair.sort()
            new = pair

        return new

    def set_description(self, set_cohps, set_coordinations_cations, set_inequivalent_cations,
                        set_infos_cation_anion_bonds, set_labels_cohps, structure):

        relevant_cations = ', '.join([str(site.specie) + str(isite) for isite, site in enumerate(structure) if isite in
                                      set_inequivalent_cations])
        self.text = []
        self.text.append(
            "The compound " + str(structure.composition.reduced_formula) + " has " + str(len(set_inequivalent_cations))
            + " symmetry-independent cations with relevant cation-anion interactions: " + relevant_cations + '.')

        for ication, ce, cation_anion_infos, labels, cohps in zip(set_inequivalent_cations, set_coordinations_cations,
                                                                  set_infos_cation_anion_bonds, set_labels_cohps,
                                                                  set_cohps):
            namecation = str(structure[ication].specie)

            # build something to get mean strengths for each type of bond
            mean_icohps = self._get_strenghts_for_each_bond(cation_anion_infos[4], cation_anion_infos[1], namecation)

            antbdg = self.get_antibdg_states(cohps, labels, namecation)
            bonds = self._get_plot_label(mean_icohps, antbdg)

            # TODO write a method that can return a dict including name of bond and True/False

            self.text.append(namecation + str(ication) + " has a " + str(self._coordination_environment_to_text(ce)
                                                                         ) + " coordination environment. It has " + str(
                bonds) + ' bonds.')

    def get_antibdg_states(self, cohps, labels, namecation):
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
    def _coordination_environment_to_text(ce):
        # TODO: include all options
        if ce == "S:1":
            return "single (CN=1)"
        if ce == "L:2":
            return "linear (CN=2)"
        if ce == "A:2":
            return "angular (CN=2)"
        if ce == "TL:3":
            return "trigonal planar (CN=3)"
        if ce == "TY:3":
            return "triangular non-coplanar (CN=3)"
        if ce == "TS:3":
            return "t-shaped (CN=3)"
        if ce == "T:4":
            return "tetrahedral (CN=4)"
        if ce == "S:4":
            return "square planar (CN=4)"
        if ce == "SY:4":
            return "square non-coplanar (CN=4)"
        if ce == "SS:4":
            return "see-saw like (CN=4)"
        if ce == "PP:5":
            return "pentagonal (CN=5)"
        if ce == "S:5":
            return "square pyramidal (CN=5)"
        if ce == "T:5":
            return "trigonal bipyramidal (CN=5)"
        if ce == "O:6":
            return "octahedral (CN=6)"
        if ce == "T:6":
            return "trigonal prismatic (CN=6)"
        if ce == "PP:6":
            return "pentagonal pyramidal (CN=6)"
        if ce == "PB:7":
            return "pentagonal bipyramidal (CN=6)"
        if ce == "ST:7":
            return "square-face capped trigonal prismatic (CN=7)"
        if ce == "ET:7":
            return "end-trigonal-face capped trigonal prismatic (CN=7)"
        if ce == "FO:7":
            return "face-capped octahedron (CN=7)"
        if ce == "C:8":
            return "cubic (CN=8)"
        if ce == "SA:8":
            return "sqaure antiprismatic (CN=8)"
        if ce == "SBT:8":
            return "square-face bicapped trigonal prismatic (CN=8)"
        if ce == "TBT:8":
            return "triangular-face bicapped trigonal prismatic (CN=8)"
        if ce == "DD:8":
            return "dodecahedronal (with triangular faces) (CN=8)"
        if ce == "DDPN":
            return "dodecahedronal (with triangular faces - p2345 plane normalized) (CN=8)"
        if ce == "HB:8":
            return "hexagonal bipyramidal (CN=8)"
        if ce == "BO_1:8":
            return "bicapped octahedral (opposed cap faces) (CN=8)"
        if ce == "BO_2:8":
            return "bicapped octahedral (cap faces with one atom in common) (CN=8)"
        if ce == "BO_3:8":
            return "bicapped octahedral (cap faces with one edge in common) (CN=8)"
        if ce == "TC:9":
            return "triangular cupola (CN=9)"
        if ce == "TT_1:9":
            return "Tricapped triangular prismatic (three square - face caps) (CN=9)"
        if ce == "TT_2:9":
            return "Tricapped triangular prismatic (two square - face caps and one triangular - face cap) (CN=9)"
        if ce == "TT_3:9":
            return "Tricapped triangular prism (one square - face cap and two triangular - face caps) (CN=9)"
        if ce == "HD:9":
            return "Heptagonal dipyramidal (CN=9)"
        if ce == "TI:9":
            return "tridiminished icosohedral (CN=9)"
        if ce == "SMA:9":
            return "Square-face monocapped antiprism (CN=9)"
        if ce == "SS:9":
            return "Square-face capped square prismatic (CN=9)"
        if ce == "TO_1:9":
            return "Tricapped octahedral (all 3 cap faces share one atom) (CN=9)"
        if ce == "TO_2:9":
            return "Tricapped octahedral (cap faces are aligned) (CN=9)"
        if ce == "TO_3:9":
            return "Tricapped octahedron (all 3 cap faces are sharing one edge of a face) (CN=9)"
        if ce == "PP:10":
            return "Pentagonal prismatic (CN=10)"
        if ce == "PA:10":
            return "Pentagonal antiprismatic (CN=10)"
        if ce == "SBSA:10":
            return "Square-face bicapped square antiprismatic (CN=10)"
        if ce == "MI:10":
            return "Metabidiminished icosahedral (CN=10)"
        if ce == "S:10":
            return "sphenocoronal (CN=10)"
        if ce == "H:10":
            return "Hexadecahedral (CN=10)"
        if ce == "BS_1:10":
            return "Bicapped square prismatic (opposite faces) (CN=10)"
        if ce == "BS_1:10":
            return "Bicapped square prismatic (opposite faces) (CN=10)"
        if ce == "BS_2:10":
            return "Bicapped square prism(adjacent faces) (CN=10)"
        if ce == "TBSA:10":
            return "Trigonal-face bicapped square antiprismatic (CN=10)"
        if ce == "PCPA:11":
            return "Pentagonal - face capped pentagonal antiprismatic (CN=11)"
        if ce == "H:11":
            return "Hendecahedral (CN=11)"
        if ce == "SH:11":
            return "Sphenoid hendecahedral (CN=11)"
        if ce == "CO:11":
            return "Cs - octahedral (CN=11)"
        if ce == "DI:11":
            return "Diminished icosahedral (CN=12)"
        if ce == "I:12":
            return "Icosahedral (CN=12)"
        if ce == "PBP: 12":
            return "Pentagonal - face bicapped pentagonal prismatic (CN=12)"
        if ce == "TT:12":
            return "Truncated tetrahedral (CN=12)"
        if ce == "C:12":
            return "Cuboctahedral (CN=12)"
        if ce == "AC:12":
            return "Anticuboctahedral (CN=12)"
        if ce == "SC:12":
            return "Square cupola (CN=12)"
        if ce == "S:12":
            return "Sphenomegacorona (CN=12)"
        if ce == "HP:12":
            return "Hexagonal prismatic (CN=12)"
        if ce == "HA:12":
            return "Hexagonal antiprismatic (CN=12)"
        if ce == "SH:13":
            return "Square-face capped hexagonal prismatic (CN=13)"

        return ce

    def get_information_all_bonds(self, chemv, lse, path_to_cohpcar, set_equivalent_sites):

        set_inequivalent_cations = []
        set_coordinations_cations = []
        set_infos_cation_anion_bonds = []
        set_labels_cohps = []
        set_cohps = []
        # only_bonds_to
        # Write something to get all different kinds of anions?
        anion_types = chemv.get_anion_types()
        for ice, ce in enumerate(lse.coordination_environments):

            # only look at inequivalent sites (use of symmetry!)!
            if ice in set_equivalent_sites and ce[0]['ce_symbol'] is not None:
                set_inequivalent_cations.append(ice)
                ce = ce[0]['ce_symbol']
                set_coordinations_cations.append(ce)
                cation_anion_infos = chemv.get_info_icohps_to_neighbors([ice])
                set_infos_cation_anion_bonds.append(cation_anion_infos)

                aniontype_labels = []
                aniontype_cohps = []

                for anion in anion_types:
                    # print(anion)
                    labels, summedcohps = chemv.get_info_cohps_to_neighbors(path_to_cohpcar, [ice],
                                                                            summed_spin_channels=True,
                                                                            per_bond=False, only_bonds_to=[str(anion)])

                    aniontype_labels.append(labels)
                    aniontype_cohps.append(summedcohps)

                set_labels_cohps.append(aniontype_labels)
                set_cohps.append(aniontype_cohps)

        return set_cohps, set_coordinations_cations, set_inequivalent_cations, set_infos_cation_anion_bonds, \
               set_labels_cohps, anion_types

    @staticmethod
    def _get_plot_label(atoms, antbd):
        plotlabels = []
        for key, item in atoms.items():
            if antbd[key] == True:
                plotlabels.append(
                    str(len(item)) + " " + str(key) + ' (mean ICOHP: ' + str(round(np.mean(item), 2)) + ' eV, '
                                                                                                        'antibonding '
                                                                                                        'interactions below EFermi)')
            else:
                plotlabels.append(
                    str(len(item)) + " " + str(key) + ' (mean ICOHP: ' + str(round(np.mean(item), 2)) + ' eV, '
                                                                                                        'no '
                                                                                                        'antibonding '
                                                                                                        'interactions below EFermi)')

        if len(plotlabels) > 1:
            plotlabel = ','.join(plotlabels[0:-1]) + ', and ' + plotlabels[-1]
        else:
            plotlabel = plotlabels[0]
        return plotlabel

    @staticmethod
    def get_bond_labels(atoms):
        all_labels = []
        for atomsnames in atoms:
            new = [LobsterNeighbors._split_string(atomsnames[0])[0], LobsterNeighbors._split_string(atomsnames[1])[0]]
            new.sort()
            # print(new2)
            string_here = new[0] + "-" + new[1]
            all_labels.append(string_here)
        count = collections.Counter(all_labels)
        return count

    def setup_env(self, path_to_icohplist, path_to_poscar):

        structure = Structure.from_file(path_to_poscar)
        sga = SpacegroupAnalyzer(structure=structure)
        equivalent_sites = sga.get_symmetry_dataset()['equivalent_atoms']
        set_equivalent_sites = list(set(equivalent_sites))
        # What do I need for an automated analysis:
        #
        self.chemv = LobsterNeighbors(filename_ICOHP=path_to_icohplist,
                                      structure=Structure.from_file(path_to_poscar),
                                      additional_condition=1,
                                      perc_strength_ICOHP=0.1
                                      )

        # determine cations and anions
        lse = self.chemv.get_light_structure_environment(only_cation_environments=True)
        return self.chemv, lse, set_equivalent_sites, structure

    def write_description(self):
        for textpart in self.text:
            print(textpart)


# TODO: write classes that allow to decribe the ICOHP situation in a crystal structure
# TODO: clean up!

# TODO: add a sentence which ICOHPs will be still considered

# TODO: einen Teil auslagern, um dann nur noch die Beschreibung hier zu haben?


class Description2:
    """
    will write a description text of all relevant bonds. It will analyse all coordination environments.
    Furthermore, it will analyse which kind of bonds contribute to the valence band close to the Fermi level

    """

    def __init__(self, analysis_object):

        self.analysis_object = analysis_object

        self.set_description()

    def set_description(self):
        self.condensed_bonding_analysis = self.analysis_object.get_condensed_bonding_analysis()
        #print(self.condensed_bonding_analysis)

        relevant_cations = ', '.join([str(site.specie) + str(isite) for isite, site in enumerate(
            self.analysis_object.structure) if isite in
                                      self.analysis_object.set_inequivalent_cations])
        self.text = []
        self.text.append("The compound " + str(self.condensed_bonding_analysis["formula"]) + " has "
                         + str(self.condensed_bonding_analysis["number_of_considered_cations"])
                         + " symmetry-independent cations with relevant cation-anion interactions: " +
                         str(relevant_cations) + '.')

        for key, item in self.condensed_bonding_analysis["sites"].items():

            # It has 3 Ta-N (mean ICOHP: -4.78 eV, antibonding interactions below EFermi),
            bond_info = []
            for type, properties in item['bonds'].items():
                if not properties["has_antibdg_states_below_Efermi"]:
                    bond_info.append(
                        str(properties['number_of_bonds']) + ' ' + item["cation"] + '-' + str(type) + ' (mean ICOHP: '
                                                                                                      '' + properties[
                            "ICOHP_mean"] + ' eV, no antibonding  interaction below EFermi)')
                else:
                    bond_info.append(
                        str(properties['number_of_bonds']) + ' ' + item["cation"] + '-' + str(type) + ' (mean ICOHP: '
                                                                                                      '' + properties[
                            "ICOHP_mean"] + ' eV, antibonding  interaction below EFermi)')

            if len(bond_info) > 1:
                bonds = ','.join(bond_info[0:-1]) + ', and ' + bond_info[-1]
            else:
                bonds = bond_info[0]
            if item["env"] == "O:6":
                self.text.append(str(item["cation"]) + str(key) + " has an " + str(
                    self._coordination_environment_to_text(item["env"]))
                                 + " coordination environment. It has " + str(bonds) + ' bonds.')
            else:
                self.text.append(str(item["cation"]) + str(key) + " has a " + str(
                    self._coordination_environment_to_text(item["env"]))
                                 + " coordination environment. It has " + str(bonds) + ' bonds.')

    def set_dict(self):
        # produce a dict with mean values for all kinds of different bonds in the structure
        self.condensed_bonding_analysis = self.analysis_object.get_condensed_bonding_analysis()
        print(self.analysis_object.list_equivalent_sites)
        print(self.analysis_object.set_inequivalent_cations)
        relevant_cation_ids = set([isite for isite in
                               self.analysis_object.list_equivalent_sites if
                               isite in self.analysis_object.set_inequivalent_cations])

        print("relevant IDs:")
        print(relevant_cation_ids)
        # self.text = []
        # self.text.append("The compound " + str(self.condensed_bonding_analysis["formula"]) + " has "
        #                  + str(self.condensed_bonding_analysis["number_of_considered_cations"])
        #                  + " symmetry-independent cations with relevant cation-anion interactions: " +
        #                  str(relevant_cations) + '.')
        #
        # output_dict={}
        final_dict = {}
        for key in relevant_cation_ids:
            item = self.condensed_bonding_analysis["sites"][key]
            # print(key)
            # print(item)
            for type, properties in item["bonds"].items():
                label = item["cation"] + '-' + str(type)
                if not label in final_dict:
                    final_dict[label] = {"number_of_bonds": int(properties["number_of_bonds"]),
                                         "ICOHP_sum": float(properties["ICOHP_sum"]),
                                         "has_antbdg": properties["has_antibdg_states_below_Efermi"]}
                else:
                    final_dict[label]["number_of_bonds"] += int(properties["number_of_bonds"])
                    final_dict[label]["ICOHP_sum"] += float(properties["ICOHP_sum"])
                    final_dict[label]["has_antbdg"] = final_dict[label]["has_antbdg"] or properties["has_antibdg_states_below_Efermi"]

        for key, item in final_dict.items():
            final_dict[key]["ICOHP_mean"] = item["ICOHP_sum"] / item["number_of_bonds"]
            final_dict[key]["number_of_bonds"] = item["number_of_bonds"]

        self.final_dict_bonds = final_dict

        final_dict_cations={}
        for key, item in self.condensed_bonding_analysis["sites"].items():

            if item["cation"] not in final_dict_cations:
                final_dict_cations[item["cation"]]=[item["env"]]
            else:
                final_dict_cations[item["cation"]].append(item["env"])
        # TODO: add coordination environment
        self.final_dict_cations=final_dict_cations

    def plot_cohps(self, save=False, filename=None):
        set_cohps = self.analysis_object.set_cohps
        set_inequivalent_cations = self.analysis_object.set_inequivalent_cations
        set_labels_cohps = self.analysis_object.set_labels_cohps
        structure = self.analysis_object.structure

        for ication, labels, cohps in zip(set_inequivalent_cations, set_labels_cohps,
                                          set_cohps):

            namecation = str(structure[ication].specie)

            cp = CohpPlotter()
            for label, cohp in zip(labels, cohps):
                if label is not None:
                    cp.add_cohp(namecation + str(ication) + ': ' + label, cohp)
            plot = cp.get_plot(integrated=False)
            plot.ylim([-4, 2])
            if save:
                plot.savefig(filename)
            plot.show()

    @staticmethod
    def _coordination_environment_to_text(ce):
        # TODO: include all options
        if ce == "S:1":
            return "single (CN=1)"
        if ce == "L:2":
            return "linear (CN=2)"
        if ce == "A:2":
            return "angular (CN=2)"
        if ce == "TL:3":
            return "trigonal planar (CN=3)"
        if ce == "TY:3":
            return "triangular non-coplanar (CN=3)"
        if ce == "TS:3":
            return "t-shaped (CN=3)"
        if ce == "T:4":
            return "tetrahedral (CN=4)"
        if ce == "S:4":
            return "square planar (CN=4)"
        if ce == "SY:4":
            return "square non-coplanar (CN=4)"
        if ce == "SS:4":
            return "see-saw like (CN=4)"
        if ce == "PP:5":
            return "pentagonal (CN=5)"
        if ce == "S:5":
            return "square pyramidal (CN=5)"
        if ce == "T:5":
            return "trigonal bipyramidal (CN=5)"
        if ce == "O:6":
            return "octahedral (CN=6)"
        if ce == "T:6":
            return "trigonal prismatic (CN=6)"
        if ce == "PP:6":
            return "pentagonal pyramidal (CN=6)"
        if ce == "PB:7":
            return "pentagonal bipyramidal (CN=6)"
        if ce == "ST:7":
            return "square-face capped trigonal prismatic (CN=7)"
        if ce == "ET:7":
            return "end-trigonal-face capped trigonal prismatic (CN=7)"
        if ce == "FO:7":
            return "face-capped octahedron (CN=7)"
        if ce == "C:8":
            return "cubic (CN=8)"
        if ce == "SA:8":
            return "sqaure antiprismatic (CN=8)"
        if ce == "SBT:8":
            return "square-face bicapped trigonal prismatic (CN=8)"
        if ce == "TBT:8":
            return "triangular-face bicapped trigonal prismatic (CN=8)"
        if ce == "DD:8":
            return "dodecahedronal (with triangular faces) (CN=8)"
        if ce == "DDPN":
            return "dodecahedronal (with triangular faces - p2345 plane normalized) (CN=8)"
        if ce == "HB:8":
            return "hexagonal bipyramidal (CN=8)"
        if ce == "BO_1:8":
            return "bicapped octahedral (opposed cap faces) (CN=8)"
        if ce == "BO_2:8":
            return "bicapped octahedral (cap faces with one atom in common) (CN=8)"
        if ce == "BO_3:8":
            return "bicapped octahedral (cap faces with one edge in common) (CN=8)"
        if ce == "TC:9":
            return "triangular cupola (CN=9)"
        if ce == "TT_1:9":
            return "Tricapped triangular prismatic (three square - face caps) (CN=9)"
        if ce == "TT_2:9":
            return "Tricapped triangular prismatic (two square - face caps and one triangular - face cap) (CN=9)"
        if ce == "TT_3:9":
            return "Tricapped triangular prism (one square - face cap and two triangular - face caps) (CN=9)"
        if ce == "HD:9":
            return "Heptagonal dipyramidal (CN=9)"
        if ce == "TI:9":
            return "tridiminished icosohedral (CN=9)"
        if ce == "SMA:9":
            return "Square-face monocapped antiprism (CN=9)"
        if ce == "SS:9":
            return "Square-face capped square prismatic (CN=9)"
        if ce == "TO_1:9":
            return "Tricapped octahedral (all 3 cap faces share one atom) (CN=9)"
        if ce == "TO_2:9":
            return "Tricapped octahedral (cap faces are aligned) (CN=9)"
        if ce == "TO_3:9":
            return "Tricapped octahedron (all 3 cap faces are sharing one edge of a face) (CN=9)"
        if ce == "PP:10":
            return "Pentagonal prismatic (CN=10)"
        if ce == "PA:10":
            return "Pentagonal antiprismatic (CN=10)"
        if ce == "SBSA:10":
            return "Square-face bicapped square antiprismatic (CN=10)"
        if ce == "MI:10":
            return "Metabidiminished icosahedral (CN=10)"
        if ce == "S:10":
            return "sphenocoronal (CN=10)"
        if ce == "H:10":
            return "Hexadecahedral (CN=10)"
        if ce == "BS_1:10":
            return "Bicapped square prismatic (opposite faces) (CN=10)"
        if ce == "BS_1:10":
            return "Bicapped square prismatic (opposite faces) (CN=10)"
        if ce == "BS_2:10":
            return "Bicapped square prism(adjacent faces) (CN=10)"
        if ce == "TBSA:10":
            return "Trigonal-face bicapped square antiprismatic (CN=10)"
        if ce == "PCPA:11":
            return "Pentagonal - face capped pentagonal antiprismatic (CN=11)"
        if ce == "H:11":
            return "Hendecahedral (CN=11)"
        if ce == "SH:11":
            return "Sphenoid hendecahedral (CN=11)"
        if ce == "CO:11":
            return "Cs - octahedral (CN=11)"
        if ce == "DI:11":
            return "Diminished icosahedral (CN=12)"
        if ce == "I:12":
            return "Icosahedral (CN=12)"
        if ce == "PBP: 12":
            return "Pentagonal - face bicapped pentagonal prismatic (CN=12)"
        if ce == "TT:12":
            return "Truncated tetrahedral (CN=12)"
        if ce == "C:12":
            return "Cuboctahedral (CN=12)"
        if ce == "AC:12":
            return "Anticuboctahedral (CN=12)"
        if ce == "SC:12":
            return "Square cupola (CN=12)"
        if ce == "S:12":
            return "Sphenomegacorona (CN=12)"
        if ce == "HP:12":
            return "Hexagonal prismatic (CN=12)"
        if ce == "HA:12":
            return "Hexagonal antiprismatic (CN=12)"
        if ce == "SH:13":
            return "Square-face capped hexagonal prismatic (CN=13)"

        return ce

    def write_description(self):
        for textpart in self.text:
            print(textpart)
