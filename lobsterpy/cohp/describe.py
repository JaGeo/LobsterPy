# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines classes to describe the COHPs automatically
"""

from lobsterpy.plotting import PlainCohpPlotter


class Description:
    """
    It will write a description text for all relevant bonds. It will analyse all relevant coordination environments in the system based on electronic structure theory.

    """

    def __init__(self, analysis_object):
        """

        Args:
            analysis_object: Analysis object from lobsterpy.analysis
        """

        self.analysis_object = analysis_object
        self.set_description()

    def set_description(self):
        """
        This class will set the description of the structures correctly.
        Important is that here the naming of the cations from the lobster files will be used
        This means that the numbers will start at 1

        Returns:

        """
        self.condensed_bonding_analysis = (
            self.analysis_object.condensed_bonding_analysis
        )
        if self.analysis_object.whichbonds == "cation-anion":
            relevant_cations = ", ".join(
                [
                    str(site.specie) + str(isite + 1)
                    for isite, site in enumerate(self.analysis_object.structure)
                    if isite in self.analysis_object.set_inequivalent_ions
                ]
            )
            self.text = []
            self.text.append(
                "The compound "
                + str(self.condensed_bonding_analysis["formula"])
                + " has "
                + str(self.condensed_bonding_analysis["number_of_considered_ions"])
                + " symmetry-independent cation(s) with relevant cation-anion interactions: "
                + str(relevant_cations)
                + "."
            )

            for key, item in self.condensed_bonding_analysis["sites"].items():

                # It has 3 Ta-N (mean ICOHP: -4.78 eV, antibonding interactions below EFermi),
                bond_info = []
                for type, properties in item["bonds"].items():
                    if not properties["has_antibdg_states_below_Efermi"]:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + " (mean ICOHP: "
                            ""
                            + properties["ICOHP_mean"]
                            + " eV, no antibonding interaction below EFermi)"
                        )
                    else:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + " (mean ICOHP: "
                            ""
                            + properties["ICOHP_mean"]
                            + " eV, antibonding interaction below EFermi)"
                        )

                if len(bond_info) > 1:
                    bonds = ",".join(bond_info[0:-1]) + ", and " + bond_info[-1]
                else:
                    bonds = bond_info[0]
                if item["env"] == "O:6":
                    self.text.append(
                        str(item["ion"])
                        + str(key + 1)
                        + " has an "
                        + str(self._coordination_environment_to_text(item["env"]))
                        + " coordination environment. It has "
                        + str(bonds)
                        + " bonds."
                    )
                else:
                    self.text.append(
                        str(item["ion"])
                        + str(key + 1)
                        + " has a "
                        + str(self._coordination_environment_to_text(item["env"]))
                        + " coordination environment. It has "
                        + str(bonds)
                        + " bonds."
                    )
        elif self.analysis_object.whichbonds == "all":
            relevant_ions = ", ".join(
                [
                    str(site.specie) + str(isite + 1)
                    for isite, site in enumerate(self.analysis_object.structure)
                    if isite in self.analysis_object.set_inequivalent_ions
                ]
            )
            self.text = []
            self.text.append(
                "The compound "
                + str(self.condensed_bonding_analysis["formula"])
                + " has "
                + str(self.condensed_bonding_analysis["number_of_considered_ions"])
                + " symmetry-independent atoms(s) with relevant bonds: "
                + str(relevant_ions)
                + "."
            )

            for key, item in self.condensed_bonding_analysis["sites"].items():

                # It has 3 Ta-N (mean ICOHP: -4.78 eV, antibonding interactions below EFermi),
                bond_info = []
                for type, properties in item["bonds"].items():
                    if not properties["has_antibdg_states_below_Efermi"]:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + " (mean ICOHP: "
                            ""
                            + properties["ICOHP_mean"]
                            + " eV, no antibonding interaction below EFermi)"
                        )
                    else:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + " (mean ICOHP: "
                            ""
                            + properties["ICOHP_mean"]
                            + " eV, antibonding interaction below EFermi)"
                        )

                if len(bond_info) > 1:
                    bonds = ",".join(bond_info[0:-1]) + ", and " + bond_info[-1]
                else:
                    bonds = bond_info[0]
                if item["env"] == "O:6":
                    self.text.append(
                        str(item["ion"])
                        + str(key + 1)
                        + " has an "
                        + str(self._coordination_environment_to_text(item["env"]))
                        + " coordination environment. It has "
                        + str(bonds)
                        + " bonds."
                    )
                else:
                    self.text.append(
                        str(item["ion"])
                        + str(key + 1)
                        + " has a "
                        + str(self._coordination_environment_to_text(item["env"]))
                        + " coordination environment. It has "
                        + str(bonds)
                        + " bonds."
                    )

        if "madelung_energy" in self.analysis_object.condensed_bonding_analysis:
            self.text.append(
                "The Madelung energy of this crystal structure per unit cell is: "
                + str(
                    self.analysis_object.condensed_bonding_analysis["madelung_energy"]
                )
                + " eV."
            )

    def plot_cohps(
        self,
        save=False,
        filename=None,
        ylim=[-4, 2],
        xlim=None,
        integrated=False,
        summed=True,
        title="",
    ):
        """
        Automatic plots of the most relevant COHP will be determined
        Args:
            save (bool): will save the plot to a file
            filename (str):
            ylim (list of float): energy scale that is shown in plot (eV)
            xlim(list of float): energy range for COHPs in eV
            integrated (bool): if True, integrated COHPs will be shown
            summed (bool): both spin cannels will be summed

        Returns:

        """
        set_cohps = self.analysis_object.set_cohps
        if self.analysis_object.whichbonds == "cation-anion":
            set_inequivalent_cations = self.analysis_object.set_inequivalent_ions
        elif self.analysis_object.whichbonds == "all":
            set_inequivalent_cations = self.analysis_object.set_inequivalent_ions
        set_labels_cohps = self.analysis_object.set_labels_cohps
        structure = self.analysis_object.structure

        for ication, labels, cohps in zip(
            set_inequivalent_cations, set_labels_cohps, set_cohps
        ):

            namecation = str(structure[ication].specie)

            cp = PlainCohpPlotter()
            for label, cohp in zip(labels, cohps):
                if label is not None:
                    cp.add_cohp(namecation + str(ication + 1) + ": " + label, cohp)
            plot = cp.get_plot(integrated=integrated)
            plot.ylim(ylim)
            if xlim is not None:
                plot.xlim(xlim)

        plot.title(title)
        if save:
            plot.savefig(filename)
        plot.show()

    @staticmethod
    def _coordination_environment_to_text(ce):
        """
        transfers a coordination environment str into a text description of the environment
        Args:
            ce (str): output from ChemEnv package (e.g., "O:6")

        Returns:
            text description of coordination environment
        """

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
            return "pentagonal bipyramidal (CN=7)"
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
        if ce == "DDPN:8":
            return (
                "dodecahedronal (with triangular faces - p2345 plane normalized) (CN=8)"
            )
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
        if ce == "1":
            return "1-fold"
        if ce == "2":
            return "2-fold"
        if ce == "3":
            return "3-fold"
        if ce == "4":
            return "4-fold"
        if ce == "5":
            return "5-fold"
        if ce == "6":
            return "6-fold"
        if ce == "7":
            return "7-fold"
        if ce == "8":
            return "8-fold"
        if ce == "9":
            return "9-fold"
        if ce == "10":
            return "10-fold"
        if ce == "11":
            return "11-fold"
        if ce == "12":
            return "12-fold"
        if ce == "13":
            return "13-fold"
        if ce == "14":
            return "14-fold"
        if ce == "15":
            return "15-fold"
        if ce == "16":
            return "16-fold"
        if ce == "17":
            return "17-fold"
        if ce == "18":
            return "18-fold"
        if ce == "19":
            return "19-fold"
        if ce == "20":
            return "20-fold"
        if ce == "21":
            return "21-fold"
        if ce == "22":
            return "22-fold"
        if ce == "23":
            return "23-fold"
        if ce == "24":
            return "24-fold"
        if ce == "25":
            return "25-fold"
        if ce == "26":
            return "26-fold"
        if ce == "27":
            return "27-fold"
        if ce == "28":
            return "28-fold"
        if ce == "29":
            return "29-fold"
        if ce == "30":
            return "30-fold"
        return ce

    def write_description(self):
        """
        prints the description of the COHPs to the screen

        """
        for textpart in self.text:
            print(textpart)
