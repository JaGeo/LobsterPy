# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module defines classes to describe the COHPs automatically."""
from __future__ import annotations

from pathlib import Path

from lobsterpy.plotting import InteractiveCohpPlotter, PlainCohpPlotter


class Description:
    """
    Base class that will write generate a text description for all relevant bonds.

    It analyses all relevant coordination environments in the system based on electronic structure theory.

    """

    def __init__(self, analysis_object):
        """
        Generate a text description for all relevant bonds.

        :param analysis_object: Analysis object from lobsterpy.analysis
        """
        self.analysis_object = analysis_object
        self.set_description()

    def set_description(self):
        """
        Set the descriptions of the structures using the cation names, starting with numbers at 1.

        Uses the cation names from the lobster files.

        Returns:
            None

        """
        self.condensed_bonding_analysis = self.analysis_object.condensed_bonding_analysis
        # set type of population analyzed
        type_pop = self.analysis_object._get_pop_type()
        # set units for populations
        units = " eV" if type_pop == "COHP" else ""
        if self.analysis_object.which_bonds == "cation-anion":
            relevant_cations = ", ".join(
                [
                    str(site.specie) + str(isite + 1)
                    for isite, site in enumerate(self.analysis_object.structure)
                    if isite in self.analysis_object.seq_ineq_ions
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
                orb_info = []
                for type, properties in item["bonds"].items():
                    if not properties["has_antibdg_states_below_Efermi"]:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + f" (mean I{type_pop}: "
                            ""
                            + properties[f"I{type_pop}_mean"]
                            + f"{units}, 0.0 percent antibonding interaction below EFermi)"
                        )
                        if self.analysis_object.orbital_resolved:
                            text_orbital = self._generate_orbital_resolved_analysis_text(
                                orbital_resolved_data=properties,
                                type_pop=type_pop,
                                atom_name=str(type),
                                ion=item["ion"],
                            )
                            orb_info.extend(text_orbital)
                    else:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + f" (mean I{type_pop}: "
                            ""
                            + properties[f"I{type_pop}_mean"]
                            + f"{units}, "
                            + str(round(properties["antibonding"]["perc"] * 100, 3))
                            + " percent antibonding interaction below EFermi)"
                        )
                        if self.analysis_object.orbital_resolved:
                            text_orbital = self._generate_orbital_resolved_analysis_text(
                                orbital_resolved_data=properties,
                                type_pop=type_pop,
                                atom_name=str(type),
                                ion=item["ion"],
                            )
                            orb_info.extend(text_orbital)

                bonds = ",".join(bond_info[0:-1]) + ", and " + bond_info[-1] if len(bond_info) > 1 else bond_info[0]

                if len(orb_info) > 1:
                    orb_bonds = "".join(orb_info).replace(".In", ". In")
                else:
                    orb_bonds = orb_info[0] if orb_info else ""
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
                    if orb_bonds:
                        self.text.append(orb_bonds)
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
                    if orb_bonds:
                        self.text.append(orb_bonds)

        elif self.analysis_object.which_bonds == "all":
            relevant_ions = ", ".join(
                [
                    str(site.specie) + str(isite + 1)
                    for isite, site in enumerate(self.analysis_object.structure)
                    if isite in self.analysis_object.seq_ineq_ions
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
                orb_info = []
                for type, properties in item["bonds"].items():
                    if not properties["has_antibdg_states_below_Efermi"]:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + f" (mean I{type_pop}: "
                            ""
                            + properties[f"I{type_pop}_mean"]
                            + f"{units}, 0.0 percent antibonding interaction below EFermi)"
                        )
                        if self.analysis_object.orbital_resolved:
                            text_orbital = self._generate_orbital_resolved_analysis_text(
                                orbital_resolved_data=properties,
                                type_pop=type_pop,
                                atom_name=str(type),
                                ion=item["ion"],
                            )
                            orb_info.extend(text_orbital)
                    else:
                        bond_info.append(
                            str(properties["number_of_bonds"])
                            + " "
                            + item["ion"]
                            + "-"
                            + str(type)
                            + f" (mean I{type_pop}: "
                            ""
                            + properties[f"I{type_pop}_mean"]
                            + f"{units}, "
                            + str(round(properties["antibonding"]["perc"] * 100, 3))
                            + " percent antibonding interaction below EFermi)"
                        )

                        if self.analysis_object.orbital_resolved:
                            text_orbital = self._generate_orbital_resolved_analysis_text(
                                orbital_resolved_data=properties,
                                type_pop=type_pop,
                                atom_name=str(type),
                                ion=item["ion"],
                            )
                            orb_info.extend(text_orbital)

                bonds = ",".join(bond_info[0:-1]) + ", and " + bond_info[-1] if len(bond_info) > 1 else bond_info[0]

                if len(orb_info) > 1:
                    orb_bonds = "".join(orb_info).replace(".In", ". In")
                else:
                    orb_bonds = orb_info[0] if orb_info else ""
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
                    if orb_bonds:
                        self.text.append(orb_bonds)
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
                    if orb_bonds:
                        self.text.append(orb_bonds)

        if (
            "madelung_energy" in self.analysis_object.condensed_bonding_analysis
            and self.analysis_object.condensed_bonding_analysis["madelung_energy"] is not None
        ):
            self.text.append(
                "The Madelung energy of this crystal structure per unit cell is: "
                + str(self.analysis_object.condensed_bonding_analysis["madelung_energy"])
                + " eV."
            )

    def _generate_orbital_resolved_analysis_text(
        self,
        orbital_resolved_data: dict,
        ion: str,
        atom_name: str,
        type_pop: str,
    ):
        """
        Generate text from orbital-resolved analysis data of the most relevant COHP, COOP, or COBI.

        :param orbital_resolved_data: dict of orbital data from condensed bonding analysis object
        :param ion: name of ion at the site
        :param atom_name: name of atomic speice to which ion is bonded
        :param type_pop: population type analysed could be "COHP" or "COOP" or "COBI"

        Returns:
            A python list with text describing the orbital which contributes
            the most to the bonding and antibonding in the bond at site
        """
        orb_info = []
        if orbital_resolved_data["orbital_data"]["orbital_summary_stats"]:
            orb_names = []
            orb_contri = []
            # get atom-pair list with ion placed first
            atom_pair = self.analysis_object._sort_name([ion, atom_name], nameion=ion)
            if "max_bonding_contribution" in orbital_resolved_data["orbital_data"]["orbital_summary_stats"]:
                for orb, data in orbital_resolved_data["orbital_data"]["orbital_summary_stats"][
                    "max_bonding_contribution"
                ].items():
                    atom_pair_with_orb_name = self.analysis_object._sort_orbital_atom_pair(
                        atom_pair=atom_pair,
                        complete_cohp=self.analysis_object.chemenv.completecohp,
                        label=orbital_resolved_data["orbital_data"]["relevant_bonds"][0],
                        orb_pair=orb,
                    )
                    orb_names.append("-".join(atom_pair_with_orb_name))
                    orb_contri.append(
                        str(
                            round(
                                data * 100,
                                3,
                            )
                        )
                    )
            orb_names_anti = []
            orb_antibonding = []
            if "max_antibonding_contribution" in orbital_resolved_data["orbital_data"]["orbital_summary_stats"]:
                for orb, data in orbital_resolved_data["orbital_data"]["orbital_summary_stats"][
                    "max_antibonding_contribution"
                ].items():
                    atom_pair_with_orb_name = self.analysis_object._sort_orbital_atom_pair(
                        atom_pair=atom_pair,
                        complete_cohp=self.analysis_object.chemenv.completecohp,
                        label=orbital_resolved_data["orbital_data"]["relevant_bonds"][0],
                        orb_pair=orb,
                    )
                    orb_names_anti.append("-".join(atom_pair_with_orb_name))
                    orb_antibonding.append(
                        str(
                            round(
                                data * 100,
                                3,
                            )
                        )
                    )
            if len(orb_contri) > 1:
                orb_name_contri = ""
                for inx, name in enumerate(orb_names):
                    if len(orb_contri) == 2 and inx + 1 != len(orb_contri):
                        orb_name_contri += f"{name} "
                    elif 2 < len(orb_contri) != inx + 1:
                        orb_name_contri += f"{name}, "
                    else:
                        orb_name_contri += f"and {name}"

                orb_name_contri += " orbitals, contributing "
                for inx, contribution in enumerate(orb_contri):
                    if len(orb_contri) == 2 and inx + 1 != len(orb_contri):
                        orb_name_contri += f"{contribution} "
                    elif 2 < len(orb_contri) != inx + 1:
                        orb_name_contri += f"{contribution}, "
                    else:
                        orb_name_contri += f"and {contribution} percent, respectively"
                num_bonds = len(orbital_resolved_data["orbital_data"]["relevant_bonds"])
                bonds = "bonds" if num_bonds > 1 else "bond"
                orb_info.append(
                    f"In the {num_bonds} "
                    + "-".join(atom_pair)
                    + f" {bonds}, relative to the summed I{type_pop}s, "
                    + "the maximum bonding contribution is from "
                    + orb_name_contri
                )
            elif not orb_contri:
                num_bonds = len(orbital_resolved_data["orbital_data"]["relevant_bonds"])
                bonds = "bonds" if num_bonds > 1 else "bond"
                orb_info.append(
                    f"In the {num_bonds} "
                    + "-".join(atom_pair)
                    + f" {bonds}, relative to the summed I{type_pop}s, "
                    + f"no orbital has a bonding contribution greater than "
                    f"{self.analysis_object.orbital_cutoff*100} percent"
                )
            else:
                num_bonds = len(orbital_resolved_data["orbital_data"]["relevant_bonds"])
                bonds = "bonds" if num_bonds > 1 else "bond"
                orb_info.append(
                    f"In the {num_bonds} "
                    + "-".join(atom_pair)
                    + f" {bonds}, relative to the summed I{type_pop}s, "
                    + "the maximum bonding contribution is from the "
                    + f"{orb_names[0]}"
                    + f" orbital, contributing {orb_contri[0]} percent"
                )

            if len(orb_antibonding) > 1:
                orb_anti = ""
                for inx, name in enumerate(orb_names_anti):
                    if len(orb_names_anti) == 2 and inx + 1 != len(orb_names_anti):
                        orb_anti += f"{name} "
                    elif 2 < len(orb_antibonding) != inx + 1:
                        orb_anti += f"{name}, "
                    else:
                        orb_anti += f"and {name}"

                orb_anti += " orbitals, contributing "
                for inx, contribution in enumerate(orb_antibonding):
                    if len(orb_names_anti) == 2 and inx + 1 != len(orb_names_anti):
                        orb_anti += f"{contribution} "
                    elif 2 < len(orb_antibonding) != inx + 1:
                        orb_anti += f"{contribution}, "
                    else:
                        orb_anti += f"and {contribution} percent, respectively."
                orb_info.append(f", whereas the maximum antibonding contribution is from {orb_anti}")
            elif not orb_antibonding:
                orb_info.append(", whereas no significant antibonding contribution is found in this bond.")
            else:
                orb_info.append(
                    f", whereas the maximum antibonding contribution is from the "
                    f"{orb_names_anti[0]} orbital, contributing {orb_antibonding[0]} percent."
                )
        else:
            # get atom-pair list with ion placed first
            atom_pair = self.analysis_object._sort_name([ion, atom_name], nameion=ion)
            percentage_cutoff = round(self.analysis_object.orbital_cutoff * 100, 2)
            orb_info.append(
                f"No individual orbital interactions detected above {percentage_cutoff} percent"
                f" with summed I{type_pop} as reference for the " + "-".join(atom_pair) + " bond."
            )

        return orb_info

    def plot_cohps(
        self,
        xlim: list[float] | None = None,
        ylim: list[float] | None = [-4, 2],
        integrated: bool = False,
        title: str = "",
        save: bool = False,
        filename: str | None = None,
        sigma: float | None = None,
        hide: bool = False,
    ):
        """
        Automatically generate plots of the most relevant COHPs, COOPs, or COBIs.

        :param save: will save the plot to a file
        :param filename: name of the file to save the plot.
        :param ylim: energy scale that is shown in plot (eV)
        :param xlim: energy range for COHPs in eV
        :param integrated: if True, integrated COHPs will be shown
        :param sigma: Standard deviation of Gaussian broadening applied to
            population data. If None, no broadening will be added.
        :param title: sets the title of figure generated
        :param hide: if True, the plot will not be shown.

        Returns:
            A matplotlib object.

        """
        seq_cohps = self.analysis_object.seq_cohps
        if self.analysis_object.which_bonds == "cation-anion":
            seq_ineq_cations = self.analysis_object.seq_ineq_ions
        elif self.analysis_object.which_bonds == "all":
            seq_ineq_cations = self.analysis_object.seq_ineq_ions
        seq_labels = self.analysis_object.seq_labels_cohps
        structure = self.analysis_object.structure

        for iplot, (ication, labels, cohps) in enumerate(zip(seq_ineq_cations, seq_labels, seq_cohps)):
            namecation = str(structure[ication].specie)

            cp = PlainCohpPlotter(
                are_coops=self.analysis_object.are_coops,
                are_cobis=self.analysis_object.are_cobis,
            )
            for label, cohp in zip(labels, cohps):
                if label is not None:
                    cp.add_cohp(namecation + str(ication + 1) + ": " + label, cohp)

            plot = cp.get_plot(integrated=integrated, sigma=sigma)
            plot.ylim(ylim)
            if xlim is not None:
                plot.xlim(xlim)

            plot.title(title)
            if save:
                if len(seq_ineq_cations) > 1:
                    if isinstance(filename, str):
                        filename = Path(filename)  # type: ignore
                    filename_new = (
                        filename.parent / f"{filename.stem}-{iplot}{filename.suffix}"  # type: ignore
                    )
                else:
                    filename_new = filename
                plot.savefig(filename_new)
        if not hide:
            plot.show()
        else:
            if not save:
                plot.close()

    def plot_interactive_cohps(
        self,
        ylim: list[float] | None = None,
        xlim: list[float] | None = None,
        save_as_html: bool = False,
        filename: str | None = None,
        integrated: bool = False,
        title: str = "",
        sigma: float | None = None,
        label_resolved: bool = False,
        orbital_resolved: bool = False,
        hide: bool = False,
    ):
        """
        Automatically generate interactive plots of the most relevant COHPs, COBIs or COOPs.

        :param save_as_html: will save the plot to a html file
        :param filename: name of the file to save the plot.
        :param ylim: energy scale that is shown in plot (eV)
        :param xlim: energy range for COHPs in eV
        :param integrated: if True, integrated COHPs will be shown
        :param sigma: Standard deviation of Gaussian broadening applied to
            population data. If None, no broadening will be added.
        :param label_resolved: if true, relevant cohp curves will be further resolved based on band labels
        :param orbital_resolved: if true, relevant orbital interactions in cohp curves will be added to figure
        :param title: Title of the interactive plot
        :param hide: if True, the plot will not be shown.

        Returns:
            A plotly.graph_objects.Figure object.
        """
        cba_cohp_plot_data = {}  # Initialize dict to store plot data
        set_cohps = self.analysis_object.seq_cohps
        set_labels_cohps = self.analysis_object.seq_labels_cohps
        set_inequivalent_cations = self.analysis_object.seq_ineq_ions
        structure = self.analysis_object.structure

        for _iplot, (ication, labels, cohps) in enumerate(zip(set_inequivalent_cations, set_labels_cohps, set_cohps)):
            label_str = f"{structure[ication].specie!s}{ication + 1!s}: "
            for label, cohp in zip(labels, cohps):
                if label is not None:
                    cba_cohp_plot_data[label_str + label] = cohp

        ip = InteractiveCohpPlotter(
            are_coops=self.analysis_object.are_coops,
            are_cobis=self.analysis_object.are_cobis,
        )
        if label_resolved or orbital_resolved:
            ip.add_all_relevant_cohps(
                analyse=self.analysis_object,
                label_resolved=label_resolved,
                orbital_resolved=orbital_resolved,
            )
        else:
            ip.add_cohps_from_plot_data(plot_data_dict=cba_cohp_plot_data)

        plot = ip.get_plot(integrated=integrated, xlim=xlim, ylim=ylim, sigma=sigma)

        plot.update_layout(title_text=title)
        if save_as_html:
            plot.write_html(filename, include_mathjax="cdn")
        if not hide:
            return plot.show()

        return plot

    @staticmethod
    def _coordination_environment_to_text(ce: str):
        """
        Convert a coordination environment string into a text description of the environment.

        :param ce: output from ChemEnv package (e.g., "O:6")

        Returns:
            A text description of coordination environment
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
            return "square antiprismatic (CN=8)"
        if ce == "SBT:8":
            return "square-face bicapped trigonal prismatic (CN=8)"
        if ce == "TBT:8":
            return "triangular-face bicapped trigonal prismatic (CN=8)"
        if ce == "DD:8":
            return "dodecahedronal (with triangular faces) (CN=8)"
        if ce == "DDPN:8":
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
        """Print the description of the COHPs or COBIs or COOPs to the screen."""
        for textpart in self.text:
            print(textpart)

    @staticmethod
    def get_calc_quality_description(quality_dict):
        """
        Generate a text description of the LOBSTER calculation quality.

        :param quality_dict: python dictionary from lobsterpy.analysis.get_lobster_calc_quality_summary
        """
        text_des = []

        for key, val in quality_dict.items():
            if key == "minimal_basis":
                if val:
                    text_des.append("The LOBSTER calculation used minimal basis.")
                if not val:
                    text_des.append(
                        "Consider rerunning the calculation with the minimum basis as well. Choosing a "
                        "larger basis set is only recommended if you see a significant improvement of "
                        "the charge spilling."
                    )

            elif key == "charge_spilling":
                text_des.append(
                    "The absolute and total charge spilling for the calculation is {} and {} %, "
                    "respectively.".format(
                        quality_dict[key]["abs_charge_spilling"],
                        quality_dict[key]["abs_total_spilling"],
                    )
                )
            elif key == "band_overlaps_analysis":
                if quality_dict[key]["file_exists"]:
                    if quality_dict[key]["has_good_quality_maxDeviation"]:
                        text_des.append(
                            "The bandOverlaps.lobster file is generated during the LOBSTER run. This "
                            "indicates that the projected wave function is not completely orthonormalized; "
                            "however, the maximal deviation values observed compared to the identity matrix "
                            "is below the threshold of 0.1."
                        )
                    else:
                        text_des.append(
                            "The bandOverlaps.lobster file is generated during the LOBSTER run. This "
                            "indicates that the projected wave function is not completely orthonormalized. "
                            "The maximal deviation value from the identity matrix is {}, and there are "
                            "{} percent k-points above the deviation threshold of 0.1. Please check the "
                            "results of other quality checks like dos comparisons, charges, "
                            "charge spillings before using the results for further "
                            "analysis.".format(
                                quality_dict[key]["max_deviation"],
                                quality_dict[key]["percent_kpoints_abv_limit"],
                            )
                        )
                else:
                    text_des.append(
                        "The projected wave function is completely orthonormalized as no "
                        "bandOverlaps.lobster file is generated during the LOBSTER run."
                    )

            elif key == "charge_comparisons":
                if val:
                    for charge in ["mulliken", "loewdin"]:
                        if val[f"bva_{charge}_agree"]:
                            text_des.append(
                                f"The atomic charge signs from {charge.capitalize()} population analysis "
                                f"agree with the bond valence analysis."
                            )
                        if not val[f"bva_{charge}_agree"]:
                            text_des.append(
                                f"The atomic charge signs from {charge.capitalize()} population analysis "
                                f"do not agree with the bond valence analysis."
                            )
                else:
                    text_des.append(
                        "Oxidation states from BVA analyzer cannot be determined. "
                        "Thus BVA charge comparison is not conducted."
                    )

            elif key == "dos_comparisons":
                comp_types = []
                tani_index = []
                for orb in val:
                    if orb.split("_")[-1] in ["s", "p", "d", "f", "summed"]:
                        comp_types.append(orb.split("_")[-1])
                        tani_index.append(str(val[orb]))
                text_des.append(
                    "The Tanimoto index from DOS comparisons in the energy range between {}, {} eV "
                    "for {} orbitals are: {}.".format(
                        val["e_range"][0],
                        val["e_range"][1],
                        ", ".join(comp_types),
                        ", ".join(tani_index),
                    )
                )

        return text_des

    @staticmethod
    def write_calc_quality_description(calc_quality_text):
        """Print the calculation quality description to the screen."""
        print(" ".join(calc_quality_text))
