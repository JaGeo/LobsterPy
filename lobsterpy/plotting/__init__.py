# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""Here classes and functions to plot Lobster outputs are provided."""

from __future__ import annotations

import typing
from itertools import cycle
from typing import TYPE_CHECKING, Any

import matplotlib as mpl
import numpy as np
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from pymatgen.core import Structure
from pymatgen.electronic_structure.cohp import Cohp, CompleteCohp, IcohpCollection
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import LobsterCompleteDos
from pymatgen.electronic_structure.plotter import CohpPlotter, DosPlotter

if TYPE_CHECKING:
    from lobsterpy.cohp.analyze import Analysis
from lobsterpy.plotting import layout_dicts as ld

base_style = resource_filename("lobsterpy.plotting", "lobsterpy_base.mplstyle")


def get_style_list(
    no_base_style: bool = False,
    styles: list[str | dict[str, Any]] | None = None,
    **kwargs,
) -> list[str | dict[str, Any]]:
    """
    Get args for matplotlib.style from user input.

    Remaining kwargs are collected as a dict and take the highest priority.

    :param no_base_style: If true, do not include lobsterpy_base.mplstyle
    :param styles: User-requested styles. These can be paths to mplstyle files,
        the names of known (matplotlib-supplied) styles, or dicts of rcParam options.
    :param kwargs: matplotlib-style sheet keyword arguments

    """
    base = [] if no_base_style else [base_style]
    if styles is None:
        styles = []

    return base + styles + [kwargs]


class PlainCohpPlotter(CohpPlotter):
    """
    Modified Pymatgen CohpPlotter with styling removed.

    This allows the styling to be manipulated more easily using matplotlib
    style sheets.

    :param zero_at_efermi: Shift all populations to have zero
            energy at the Fermi level. Defaults to True.
    :param are_coops: Bool indicating that populations are COOPs, not COHPs.
                Defaults to False for COHPs.
    :param are_cobis: Bool indicating that populations are COBIs, not COHPs.
                Defaults to False for COHPs.
    """

    def get_plot(
        self,
        ax: mpl.axes.Axes | None = None,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        plot_negative: bool | None = None,
        integrated: bool = False,
        invert_axes: bool = True,
        sigma: float | None = None,
    ):
        """
        Get a matplotlib plot showing the COHP.

        :param ax: Existing Matplotlib Axes object to plot to.
        :param xlim: Specifies the x-axis limits. Defaults to None for
            automatic determination.
        :param ylim: Specifies the y-axis limits. Defaults to None for
            automatic determination.
        :param plot_negative: It is common to plot -COHP(E) so that the
            sign means the same for COOPs and COHPs. Defaults to None
            for automatic determination: If are_coops is True, this
            will be set to False, else it will be set to True.
        :param integrated: Switch to plot ICOHPs. Defaults to False.
        :param invert_axes: Put the energies onto the y-axis, which is
            common in chemistry.
        :param sigma: Standard deviation of Gaussian broadening applied to
            population data. If this is unset (None) no broadening will be
            added.

        Returns:
            A matplotlib object.

        """
        if self.are_coops and not self.are_cobis:
            cohp_label = "COOP"
        elif self.are_cobis and not self.are_coops:
            cohp_label = "COBI"
        elif self.are_cobis and self.are_coops:
            raise ValueError("Plot data should not contain COBI and COOP data at same time")
        else:
            cohp_label = "COHP" + " (eV)"

        if plot_negative is None:
            plot_negative = (not self.are_coops) and (not self.are_cobis)

        if integrated:
            cohp_label = "I" + cohp_label

        if plot_negative:
            cohp_label = "$-$" + cohp_label

        energy_label = "$E - E_f$ (eV)" if self.zero_at_efermi else "$E$ (eV)"

        colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
        ncolors = len(colors)

        if ax is None:
            _, ax = plt.subplots()

        allpts = []
        keys = self._cohps.keys()
        for i, key in enumerate(keys):
            energies = self._cohps[key]["energies"]
            populations = self._cohps[key]["COHP"] if not integrated else self._cohps[key]["ICOHP"]
            for spin in [Spin.up, Spin.down]:
                if spin in populations:
                    if invert_axes:
                        x = -populations[spin] if plot_negative else populations[spin]
                        y = energies
                        x = self._broaden(y, x, sigma=sigma)
                    else:
                        x = energies
                        y = -populations[spin] if plot_negative else populations[spin]
                        y = self._broaden(x, y, sigma=sigma)
                    allpts.extend(list(zip(x, y)))
                    if spin == Spin.up:
                        ax.plot(
                            x,
                            y,
                            color=colors[i % ncolors],
                            linestyle="-",
                            label=str(key),
                        )
                    else:
                        ax.plot(x, y, color=colors[i % ncolors], linestyle="--")

        if xlim:
            ax.set_xlim(xlim)
        xlim = ax.get_xlim()
        assert isinstance(xlim, tuple)

        if ylim:
            ax.set_ylim(ylim)
        else:
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            if relevanty:
                plt.ylim((min(relevanty), max(relevanty)))

        grid_like_line_kwargs = {
            "color": mpl.rcParams["grid.color"],
            "linewidth": mpl.rcParams["grid.linewidth"],
            "linestyle": mpl.rcParams["grid.linestyle"],
            "alpha": mpl.rcParams["grid.alpha"],
            "zorder": 0,
        }

        if not invert_axes:
            ax.axhline(**grid_like_line_kwargs)

            if self.zero_at_efermi:
                ax.axvline(**grid_like_line_kwargs)

            else:
                ax.axvline(self._cohps[key]["efermi"], **grid_like_line_kwargs)
        else:
            ax.axvline(**grid_like_line_kwargs)

            if self.zero_at_efermi:
                ax.axhline(**grid_like_line_kwargs)
            else:
                ax.axhline(self._cohps[key]["efermi"], **grid_like_line_kwargs)

        if invert_axes:
            plt.xlabel(cohp_label)
            plt.ylabel(energy_label)
        else:
            plt.xlabel(energy_label)
            plt.ylabel(cohp_label)

        _ = ax.legend()

        return plt

    @staticmethod
    def _broaden(energies: np.ndarray, population: np.ndarray, sigma=None, cutoff=4.0):
        """
        Broaden the spectrum with a given standard deviation.

        The population is convolved with a normalised Gaussian kernel. This
        requires the energy grid to be regularly-spaced.

        :param energies: Regularly-spaced energy series
        :param population: Population data for broadening
        :param sigma: Standard deviation for Gaussian broadening. If sigma is None
            then the input data is returned without any processing.
        :param cutoff: Range cutoff for broadening kernel, as a multiple of sigma.

        Return:
            Broadened population

        """
        from scipy.signal import convolve
        from scipy.stats import norm

        if sigma is None:
            return population

        spacing = np.mean(np.diff(energies))
        if not np.allclose(np.diff(energies), spacing, atol=1e-5):
            raise ValueError("Energy grid is not regular, cannot broaden with discrete convolution.")

        # Obtain symmetric mesh for broadening kernel, centered on zero
        kernel_x = np.arange(0, cutoff * sigma + 0.5 * spacing, spacing)
        kernel_x = np.concatenate([-kernel_x[-1:1:-1], kernel_x])

        kernel = norm.pdf(kernel_x, scale=sigma)

        return convolve(population, kernel, mode="same") / kernel.sum()


class PlainDosPlotter(DosPlotter):
    """
    Modified Pymatgen DosPlotter with styling removed.

    This allows the styling to be manipulated more easily using matplotlib
    style sheets. It also adds additional functionalities to plotter

    :param zero_at_efermi: Shift all DOS to have zero
            energy at the Fermi level. Defaults to True.
    :param stack: Bool indicating that plot should be stacked
            area graph.
    :param sigma: Standard deviation for gaussian smearing.
    :param summed: Will plot summed dos spin populations.
            Defaults to False.
    """

    def __init__(
        self,
        zero_at_efermi: bool = True,
        stack: bool = False,
        summed: bool = False,
        sigma: float | None = None,
    ) -> None:
        """
        Initialize DOS plotter.

        :param zero_at_efermi: Whether to shift all Dos to have zero energy at the
            fermi energy. Defaults to True.
        :param stack: Whether to plot the DOS as a stacked area graph
        :param summed: Whether to plot the summed spins DOS.
        :param sigma: Specify a standard deviation for Gaussian smearing
            the DOS for nicer looking plots. Defaults to None for no smearing.

        """
        super().__init__(zero_at_efermi, stack, sigma)
        self.zero_at_efermi = zero_at_efermi
        self.stack = stack
        self.sigma = sigma
        self.summed = summed
        self._norm_val = True
        self._doses = {}  # type: ignore

    def add_dos(self, label: str, dos: LobsterCompleteDos) -> None:
        """
        Add a dos for plotting.

        :param label: label for the DOS. Must be unique.
        :param dos: LobsterCompleteDos object

        """
        if dos.norm_vol is None:
            self._norm_val = False
        energies = dos.energies
        if self.summed:
            if self.sigma:
                smeared_densities = dos.get_smeared_densities(self.sigma)
                if Spin.down in smeared_densities:
                    added_densities = smeared_densities[Spin.up] + smeared_densities[Spin.down]
                    densities = {Spin.up: added_densities}
                else:
                    densities = smeared_densities
            else:
                densities = {Spin.up: dos.get_densities()}
        else:
            densities = dos.get_smeared_densities(self.sigma) if self.sigma else dos.densities

        efermi = dos.efermi

        self._doses[label] = {
            "energies": energies,
            "densities": densities,
            "efermi": efermi,
        }

    def add_site_orbital_dos(self, dos: LobsterCompleteDos, orbital: str, site_index: int):
        """
        Add orbital dos at particular site.

        :param dos: LobsterCompleteDos object
        :param orbital: Orbitals name at the site. Must be unique.
        :param site_index: site index in the structure

        """
        if dos.norm_vol is None:
            self._norm_val = False
        site = dos.structure.sites[site_index]
        avail_orbs = list(dos.pdos[site])
        if orbital not in avail_orbs and orbital != "all":
            str_orbs = ", ".join(avail_orbs)
            raise ValueError(f"Requested orbital is not available for this site, available orbitals are {str_orbs}")

        if orbital == "all":
            for orb in avail_orbs:
                dos_obj = dos.get_site_orbital_dos(site=site, orbital=orb)
                label = site.species_string + str(site_index + 1) + f": {orb}"
                energies = dos_obj.energies
                if self.summed:
                    if self.sigma:
                        smeared_densities = dos_obj.get_smeared_densities(self.sigma)
                        if Spin.down in smeared_densities:
                            added_densities = smeared_densities[Spin.up] + smeared_densities[Spin.down]
                            densities = {Spin.up: added_densities}
                        else:
                            densities = smeared_densities
                    else:
                        densities = {Spin.up: dos_obj.get_densities()}
                else:
                    densities = dos_obj.get_smeared_densities(self.sigma) if self.sigma else dos_obj.densities

                efermi = dos_obj.efermi

                self._doses[label] = {
                    "energies": energies,
                    "densities": densities,
                    "efermi": efermi,
                }
        else:
            dos_obj = dos.get_site_orbital_dos(site=site, orbital=orbital)
            label = site.species_string + str(site_index + 1) + f": {orbital}"

            energies = dos_obj.energies
            if self.summed:
                if self.sigma:
                    smeared_densities = dos_obj.get_smeared_densities(self.sigma)
                    if Spin.down in smeared_densities:
                        added_densities = smeared_densities[Spin.up] + smeared_densities[Spin.down]
                        densities = {Spin.up: added_densities}
                    else:
                        densities = smeared_densities
                else:
                    densities = {Spin.up: dos_obj.get_densities()}
            else:
                densities = dos_obj.get_smeared_densities(self.sigma) if self.sigma else dos_obj.densities

            efermi = dos_obj.efermi

            self._doses[label] = {
                "energies": energies,
                "densities": densities,
                "efermi": efermi,
            }

    @typing.no_type_check
    def get_plot(
        self,
        ax: mpl.axes.Axes | None = None,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        invert_axes: bool = False,
        beta_dashed: bool = False,
        sigma: float | None = None,
    ):
        """
        Get a matplotlib plot showing the COHP.

        :param ax: Existing Matplotlib Axes object to plot to.
        :param xlim: Specifies the x-axis limits. Defaults to None for
            automatic determination.
        :param ylim: Specifies the y-axis limits. Defaults to None for
            automatic determination.
        :param invert_axes: Put the energies onto the y-axis, which is
            common in chemistry.
        :param beta_dashed: Plots the beta spin channel with a dashed line. Defaults to False
        :param sigma: Standard deviation of Gaussian broadening applied to population data.
            If this is unset (None) no broadening will be added.

        Returns:
            A matplotlib object.

        """
        ys = None
        all_densities = []
        all_energies = []

        colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
        n_colors = len(colors)

        if ax is None:
            _, ax = plt.subplots()

        # Note that this complicated processing of energies is to allow for
        # stacked plots in matplotlib.
        for dos in self._doses.values():
            energies = dos["energies"]
            densities = dos["densities"]
            if not ys:
                ys = {
                    Spin.up: np.zeros(energies.shape),
                    Spin.down: np.zeros(energies.shape),
                }
            new_dens = {}
            for spin in [Spin.up, Spin.down]:
                if spin in densities:
                    if self.stack:
                        ys[spin] += densities[spin]
                        new_dens[spin] = ys[spin].copy()
                    else:
                        new_dens[spin] = densities[spin]
            all_energies.append(energies)
            all_densities.append(new_dens)

        keys = list(reversed(self._doses))
        all_densities.reverse()
        all_energies.reverse()
        all_pts = []

        for idx, key in enumerate(keys):
            for spin in [Spin.up, Spin.down]:
                if spin in all_densities[idx]:
                    energy = all_energies[idx]
                    densities = list(int(spin) * all_densities[idx][spin])
                    if invert_axes:
                        x = densities
                        y = energy
                    else:
                        x = energy
                        y = densities
                    all_pts.extend(list(zip(x, y)))
                    if self.stack:
                        ax.fill(x, y, color=colors[idx % n_colors], label=str(key))
                    elif spin == Spin.down and beta_dashed:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % n_colors],
                            label=str(key),
                            linestyle="--",
                            linewidth=3,
                        )
                    else:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % n_colors],
                            label=str(key),
                            linewidth=3,
                        )

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        elif not invert_axes:
            xlim = ax.get_xlim()
            relevant_y = [p[1] for p in all_pts if xlim[0] < p[0] < xlim[1]]
            ax.set_ylim((min(relevant_y), max(relevant_y)))
        if not xlim and invert_axes:
            ylim = ax.get_ylim()
            relevant_y = [p[0] for p in all_pts if ylim[0] < p[1] < ylim[1]]
            ax.set_xlim((min(relevant_y), max(relevant_y)))

        if self.zero_at_efermi:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            (ax.plot(xlim, [0, 0], "k--", linewidth=2) if invert_axes else ax.plot([0, 0], ylim, "k--", linewidth=2))

        if invert_axes:
            ax.set_ylabel("Energies (eV)")
            ax.set_xlabel(f"Density of states (states/eV{'/Å³' if self._norm_val else ''})")
            ax.axvline(x=0, color="k", linestyle="--", linewidth=2)
        else:
            ax.set_xlabel("Energies (eV)")
            if self._norm_val:
                ax.set_ylabel("Density of states (states/eV/Å³)")
            else:
                ax.set_ylabel("Density of states (states/eV)")
            ax.axhline(y=0, color="k", linestyle="--", linewidth=2)

        # Remove duplicate labels with a dictionary
        handles, labels = ax.get_legend_handles_labels()
        label_dict = dict(zip(labels, handles))
        ax.legend(label_dict.values(), label_dict.keys())
        legend_text = ax.get_legend().get_texts()  # all the text.Text instance in the legend
        plt.setp(legend_text, fontsize=30)
        plt.tight_layout()
        _ = ax.legend()

        return plt


class InteractiveCohpPlotter(CohpPlotter):
    """Interactive COHP, COBI or COOP plotter to view all relevant bonds in one figure.

    :param zero_at_efermi: Shift all populations to have zero
            energy at the Fermi level. Defaults to True.
    :param are_coops: Bool indicating that populations are COOPs, not COHPs.
                Defaults to False for COHPs.
    :param are_cobis: Bool indicating that populations are COBIs, not COHPs.
                Defaults to False for COHPs.
    """

    COLOR_PALETTE = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#ffff33",
        "#a65628",
        "#f781bf",
        "#999999",
    ]

    def add_cohp(self, label: str, cohp: Cohp):
        """
        Add COHP object to the plotter.

        :param label: Label for the COHP. Must be unique.
        :param cohp: COHP object.

        """
        if "All" not in self._cohps:
            self._cohps["All"] = {}

        energies = cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies

        if label not in self._cohps["All"]:
            self._cohps["All"].update(
                {
                    label: {
                        "energies": energies,
                        "COHP": cohp.get_cohp(),
                        "ICOHP": cohp.get_icohp(),
                        "efermi": cohp.efermi,
                    }
                }
            )
        else:
            raise ValueError(
                "Please use another label to add the COHP, provided label already exists "
                "in the plot data, which will result in overwriting the existing COHP data."
            )

    def add_all_relevant_cohps(
        self,
        analyse: Analysis,
        suffix: str = "",
        label_resolved: bool = True,
        orbital_resolved: bool = False,
    ) -> None:
        """
        Add all relevant COHPs from lobsterpy analyse object.

        :param analyse: Analyse object from lobsterpy.
        :param suffix: Optional addition to LOBSTER label to avoid key conflicts when plotting multiple
            calcs or just for additional legend information.
        :param label_resolved: bool indicating to obtain label resolved interactive plots for relevant bonds.
            If false, will return summed cohp curves of unique relevant bonds.
        :param orbital_resolved: bool indicating to include orbital resolved interactive cohps for relevant bonds.

        """
        complete_cohp = analyse.chemenv.completecohp

        plot_data = analyse.get_site_bond_resolved_labels()

        if "All" not in self._cohps:
            self._cohps["All"] = {}

        # iterate and extract the data to be plotted from cohp objects
        for bond_key, labels in plot_data.items():
            count = len(labels)
            label_with_count = self._insert_number_of_bonds_in_label(
                label=bond_key, character=":", number_of_bonds=count
            )
            # get summed cohps from the relevant bond label at the site
            cohp = complete_cohp.get_summed_cohp_by_label_list(label_list=labels)
            energies = cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies
            drop_down_key = plot_legend = label_with_count + suffix  # label for the dropdown menu
            self._cohps["All"].update(
                {
                    plot_legend: {
                        "energies": energies,
                        "COHP": cohp.get_cohp(),
                        "ICOHP": cohp.get_icohp(),
                        "efermi": cohp.efermi,
                    }
                }
            )
            if len(plot_data) > 1:
                if drop_down_key not in self._cohps:
                    self._cohps[drop_down_key] = {}
                self._cohps[drop_down_key].update(
                    {
                        plot_legend: {
                            "energies": energies,
                            "COHP": cohp.get_cohp(),
                            "ICOHP": cohp.get_icohp(),
                            "efermi": cohp.efermi,
                        }
                    }
                )

            # Add cohp data for each relevant bond label iteratively
            if label_resolved and not orbital_resolved:
                if label_with_count + suffix not in self._cohps:
                    self._cohps[label_with_count + suffix] = {}
                # Get cohp data for each relevant bond label at the site
                for label in labels:
                    cohp = complete_cohp.get_cohp_by_label(label)
                    energies = cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies
                    drop_down_key = label_with_count + suffix  # label for the dropdown menu
                    plot_legend = self._get_plot_label_for_label_resolved(
                        structure=analyse.structure,
                        label_list=[label],
                        complete_cohp=complete_cohp,
                        orb_list=[],
                        label_resolved=True,
                        orbital_resolved=False,
                    )
                    if len(plot_data) > 1:
                        self._cohps[drop_down_key].update(
                            {
                                plot_legend: {
                                    "energies": energies,
                                    "COHP": cohp.get_cohp(),
                                    "ICOHP": cohp.get_icohp(),
                                    "efermi": cohp.efermi,
                                }
                            }
                        )

                    plot_legend_here = plot_legend + suffix
                    self._cohps["All"].update(
                        {
                            plot_legend_here: {
                                "energies": energies,
                                "COHP": cohp.get_cohp(),
                                "ICOHP": cohp.get_icohp(),
                                "efermi": cohp.efermi,
                            }
                        }
                    )
            # Adds cohp data for each relevant orbitals and bond label iteratively
            if orbital_resolved and label_resolved:
                if label_with_count + suffix not in self._cohps:
                    self._cohps[label_with_count + suffix] = {}
                # get relevant orbitals associated with each bond label
                plot_data_orb = analyse.get_site_orbital_resolved_labels()
                drop_down_key = label_with_count + suffix  # label for the dropdown menu
                key_val = plot_data_orb[bond_key]
                # get cohp data for each orbital and associated bond labels iteratively
                for orb, val in key_val.items():
                    for lab in val:
                        cohp_orb = complete_cohp.get_summed_cohp_by_label_and_orbital_list(
                            label_list=[lab], orbital_list=[orb]
                        )

                        energies = cohp_orb.energies - cohp_orb.efermi if self.zero_at_efermi else cohp_orb.energies
                        # plot legends will contain species and orbital along with relevant bond label
                        plot_legend = self._get_plot_label_for_label_resolved(
                            structure=analyse.structure,
                            label_list=[lab],
                            complete_cohp=complete_cohp,
                            orb_list=[orb],
                            orbital_resolved=True,
                            label_resolved=True,
                        )
                        if len(plot_data) > 1:
                            self._cohps[drop_down_key].update(
                                {
                                    plot_legend: {
                                        "energies": energies,
                                        "COHP": cohp_orb.get_cohp(),
                                        "ICOHP": cohp_orb.get_icohp(),
                                        "efermi": cohp_orb.efermi,
                                    }
                                }
                            )

                        plot_legend_here = plot_legend + suffix

                        self._cohps["All"].update(
                            {
                                plot_legend_here: {
                                    "energies": energies,
                                    "COHP": cohp_orb.get_cohp(),
                                    "ICOHP": cohp_orb.get_icohp(),
                                    "efermi": cohp_orb.efermi,
                                }
                            }
                        )
            # Adds summed cohp data for each relevant orbitals
            if orbital_resolved and not label_resolved:
                if label_with_count + suffix not in self._cohps:
                    self._cohps[label_with_count + suffix] = {}
                # get relevant orbitals associated with each bond label
                plot_data_orb = analyse.get_site_orbital_resolved_labels()
                drop_down_key = label_with_count + suffix  # label for the dropdown menu
                key_val = plot_data_orb[bond_key]
                # get summed cohp data for each relevant orbital
                for orb, val in key_val.items():
                    cohp_orb = complete_cohp.get_summed_cohp_by_label_and_orbital_list(
                        label_list=val, orbital_list=[orb] * len(val)
                    )

                    energies = cohp_orb.energies - cohp_orb.efermi if self.zero_at_efermi else cohp_orb.energies
                    # plot legends will contain species and orbital along with number of bonds at the site
                    plot_legend = self._get_plot_label_for_label_resolved(
                        structure=analyse.structure,
                        label_list=val,
                        complete_cohp=complete_cohp,
                        orb_list=[orb],
                        orbital_resolved=True,
                        label_resolved=False,
                    )

                    if len(plot_data) > 1:
                        self._cohps[drop_down_key].update(
                            {
                                plot_legend: {
                                    "energies": energies,
                                    "COHP": cohp_orb.get_cohp(),
                                    "ICOHP": cohp_orb.get_icohp(),
                                    "efermi": cohp_orb.efermi,
                                }
                            }
                        )

                    plot_legend_here = plot_legend + suffix

                    self._cohps["All"].update(
                        {
                            plot_legend_here: {
                                "energies": energies,
                                "COHP": cohp_orb.get_cohp(),
                                "ICOHP": cohp_orb.get_icohp(),
                                "efermi": cohp_orb.efermi,
                            }
                        }
                    )

    def add_cohps_by_lobster_label(self, analyse: Analysis, label_list: list[str], suffix: str = ""):
        """
        Add COHPs explicitly specified in label list.

        :param analyse: Analyse object from lobsterpy.
        :param label_list: List of COHP labels as from LOBSTER.
        :param suffix: Optional addition to LOBSTER label to avoid key
            conflicts when plotting multiple calcs or just for additional legend information.

        """
        complete_cohp = analyse.chemenv.completecohp

        if "All" not in self._cohps:
            self._cohps["All"] = {}

        for label in label_list:
            atom1 = complete_cohp.bonds[label]["sites"][0].species_string
            atom2 = complete_cohp.bonds[label]["sites"][1].species_string
            sorted_label = sorted([atom1, atom2])
            cohp = complete_cohp.get_cohp_by_label(label)
            energies = cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies
            key = sorted_label[0] + "-" + sorted_label[1] + ": " + label + suffix
            self._cohps["All"].update(
                {
                    key: {
                        "energies": energies,
                        "COHP": cohp.get_cohp(),
                        "ICOHP": cohp.get_icohp(),
                        "efermi": cohp.efermi,
                    }
                }
            )

    def add_cohps_from_plot_data(self, plot_data_dict: dict[str, Cohp], suffix: str = ""):
        """
        Add all relevant COHPs for specified bond type from lobster lightweight json.gz file.

        :param plot_data_dict: Lobsterpy plot data dict
        :param suffix: Optional addition to LOBSTER label to avoid key
            conflicts when plotting multiple calcs or just for additional legend information.

        """
        # convert to cohp objects
        plot_data_dict = plot_data_dict.copy()
        for key, cohps in plot_data_dict.items():
            if isinstance(cohps, Cohp):
                plot_data_dict.update({key: cohps})
            else:
                try:
                    cohps = Cohp.from_dict(cohps)
                    plot_data_dict.update({key: cohps})
                except (TypeError, AttributeError):
                    raise ValueError(
                        "The data provided could not be converted to cohp object.Please recheck the input data"
                    )

        if "All" not in self._cohps:
            self._cohps["All"] = {}

        for bond_key, cohps in plot_data_dict.items():
            energies = cohps.energies - cohps.efermi if self.zero_at_efermi else cohps.energies
            key = bond_key + suffix
            self._cohps["All"].update(
                {
                    key: {
                        "energies": energies,
                        "COHP": cohps.get_cohp(),
                        "ICOHP": cohps.get_icohp(),
                        "efermi": cohps.efermi,
                    }
                }
            )

    def get_plot(
        self,
        xlim: list[float] | None = None,
        rangeslider: bool = False,
        ylim: list[float] | None = None,
        plot_negative: bool | None = None,
        integrated: bool = False,
        invert_axes: bool = True,
        sigma: float | None = None,
        colors: list[str] | None = None,
    ):
        """
        Get an interactive plotly figure showing the COHPs.

        :param xlim: Specifies the x-axis limits. Defaults to None for
            automatic determination.
        :param rangeslider: Adds a plotly.graph_objs.layout.xaxis.Rangeslider
            object to figure to allow easy manipulation of x-axis..
        :param ylim: Specifies the y-axis limits. Defaults to None for
            automatic determination.
        :param plot_negative: It is common to plot -COHP(E) so that the
            sign means the same for COOPs and COHPs. Defaults to None
            for automatic determination - If are_coops is True, this
            will be set to False, else it will be set to True.
        :param integrated: Switch to plot ICOHPs. Defaults to False.
        :param invert_axes: Put the energies onto the y-axis, which is
            common in chemistry.
        :param sigma: Standard deviation of Gaussian broadening applied to
            population data. If this is unset (None) no broadening will be added.
        :param colors: list of hex color codes to be used in plot

        Returns:
            A  plotly.graph_objects.Figure object.

        """
        if self.are_coops and not self.are_cobis:
            cohp_label = "COOP"
        elif self.are_cobis and not self.are_coops:
            cohp_label = "COBI"
        elif self.are_cobis and self.are_coops:
            raise ValueError("Plot data should not contain COBI and COOP data at same time")
        else:
            cohp_label = "COHP" + " (eV)"

        if plot_negative is None:
            plot_negative = (not self.are_coops) and (not self.are_cobis)

        if integrated:
            cohp_label = "I" + cohp_label

        if plot_negative:
            cohp_label = "\u2212" + cohp_label

        energy_label = "$E - E_f \\text{ (eV)}$" if self.zero_at_efermi else "$E \\text{ (eV)}$"

        # Setting up repeating color scheme (same as for matplotlib plots in .mplstyle)
        palette = InteractiveCohpPlotter.COLOR_PALETTE if colors is None else colors
        pal_iter = cycle(palette)

        traces = {}  # type: ignore
        for k, v in self._cohps.items():
            traces.update({k: []})
            for label in v:
                population_key = v[label]["ICOHP"] if integrated else v[label]["COHP"]
                band_color = next(pal_iter)
                for spin in [Spin.up, Spin.down]:
                    if spin in population_key:
                        population = [-i for i in population_key[spin]] if plot_negative else population_key[spin]
                        if invert_axes:
                            x = population
                            y = v[label]["energies"]
                            x = PlainCohpPlotter._broaden(y, x, sigma=sigma)
                        else:
                            x = v[label]["energies"]
                            y = population
                            y = PlainCohpPlotter._broaden(x, y, sigma=sigma)
                        if spin == Spin.up:
                            trace = go.Scatter(x=x, y=y, name=label)
                            trace.update(ld.spin_up_trace_style_dict)
                        else:
                            trace = go.Scatter(x=x, y=y, name="")
                            trace.update(ld.spin_down_trace_style_dict)
                        trace.update(line={"color": band_color})
                        traces[k].append(trace)

        energy_axis = (
            go.layout.YAxis(title=energy_label)
            if invert_axes
            else go.layout.XAxis(title=energy_label, rangeslider={"visible": rangeslider})
        )
        energy_axis.update(ld.energy_axis_style_dict)
        cohp_axis = (
            go.layout.XAxis(title=cohp_label, rangeslider={"visible": rangeslider})
            if invert_axes
            else go.layout.YAxis(title=cohp_label)
        )
        cohp_axis.update(ld.cohp_axis_style_dict)

        layout = (
            go.Layout(xaxis=cohp_axis, yaxis=energy_axis)
            if invert_axes
            else go.Layout(xaxis=energy_axis, yaxis=cohp_axis)
        )

        # Create figure object
        fig = go.Figure(layout=layout)
        fig.update_layout(ld.layout_dict)
        fig.update_layout(legend=ld.legend_style_dict)

        # Add all traces to figure
        for val_trace in traces.values():
            for trace in val_trace:
                fig.add_trace(trace)

        # Update layout with dropdown menu if label resolved plot
        if len(traces) > 2:
            # Update visibility of traces
            for i, _ in enumerate(fig.data):
                if i <= len(traces["All"]) - 1:
                    pass
                else:
                    fig.data[i].visible = False
            # Add dropdown buttons
            fig.update_layout(
                updatemenus=[
                    {
                        "buttons": [
                            {
                                "args": [
                                    {
                                        "visible": [
                                            selected_group == group
                                            for group, val_trace in traces.items()
                                            for _trace in val_trace
                                        ]
                                    }
                                ],
                                "label": selected_group,
                                "method": "update",
                            }
                            for selected_group in traces
                        ],
                        "direction": "down",
                        "showactive": True,
                        "active": 0,
                        "x": 0.5,
                        "y": 1.15,
                        "bgcolor": "rgba(255,255,255,0.8)",
                        "bordercolor": "rgba(0,0,0,0.2)",
                        "xanchor": "center",
                        "yanchor": "top",
                        "font": {"family": "Arial", "color": "#444444", "size": 18},
                    }
                ]
            )

        if xlim:
            fig.update_xaxes(range=xlim)
        if ylim:
            fig.update_yaxes(range=ylim)

        fig.update_yaxes(automargin=True)

        return fig

    @staticmethod
    def _insert_number_of_bonds_in_label(label: str, character: str, number_of_bonds: int) -> str:
        """
        Add number of bonds to bond label.

        For example : for input label 'Ba1: Ba-Ti', character ':', number_of_bonds: 3,
        Will return 'Ba1: 3 x Ba-Ti'

        :param label: bond label to which number of bonds needs to be inserted
        :param character: string character where number of bonds needs to be inserted
        :param number_of_bonds: number of bonds corresponding to the label

        Returns:
             bond label with number of bonds inserted
        """
        return label.replace(character, f"{character} {number_of_bonds} x", 1)

    @staticmethod
    def _get_plot_label_for_label_resolved(
        structure: Structure,
        label_list: list[str],
        complete_cohp: CompleteCohp,
        orb_list: list[str],
        label_resolved: bool = False,
        orbital_resolved: bool = False,
    ) -> str:
        """
        Generate plot labels for both orbital and label-resolved plots.

        For example for NaCl structure, label:21, orb: 3s-3s
        Will return '21: Cl2-Na1 (2.85 Å)' if label_resolved is True and orbital_resolved is False
        Will return '21: Cl2(3s)-Na1(3s) (2.85 Å)' If label and orbital resolved True
        Will return 'Cl(3s)-Na(3s) (2.85 Å)' if orbital_resolved is True and label_resolved is False

        :param structure: pymatgen structure object
        :param label_list: bond label to which number of bonds needs to be inserted
        :param complete_cohp:  complete cohp object
        :param orb_list: relevant orbital
        :param label_resolved: specifies type of label to be returned is for label_resolved case
        :param orbital_resolved: specifies type of label to be returned is for orbital_resolved case

        Returns:
             plot label string
        """
        if label_resolved and not orbital_resolved:
            atom_pairs = []
            for site in complete_cohp.bonds[label_list[0]]["sites"]:
                atom = f"{site.species_string}{structure.sites.index(site) + 1!s}"
                atom_pairs.append(atom)
            atom_pair_str = "-".join(atom_pairs)
            bond_length = round(complete_cohp.bonds[label_list[0]]["length"], 2)
            plot_label = f"{label_list[0]}: {atom_pair_str} ({bond_length} \u00c5)"

        elif not label_resolved and orbital_resolved:
            orb_atom_pairs = []
            orb_pair = orb_list[0].split("-")
            for site, site_orb in zip(complete_cohp.bonds[label_list[0]]["sites"], orb_pair):
                atom_orb = f"{site.species_string} ({site_orb})"
                orb_atom_pairs.append(atom_orb)
            orb_atom_pairs_str = "-".join(orb_atom_pairs)
            bond_length = round(complete_cohp.bonds[label_list[0]]["length"], 2)
            plot_label = f"{len(label_list)}x {orb_atom_pairs_str} ({bond_length} \u00c5)"
        else:
            orb_atom_pairs = []
            orb_pair = orb_list[0].split("-")
            for site, site_orb in zip(complete_cohp.bonds[label_list[0]]["sites"], orb_pair):
                atom_orb = f"{site.species_string}{structure.sites.index(site) + 1!s} ({site_orb})"
                orb_atom_pairs.append(atom_orb)
            orb_atom_pairs_str = "-".join(orb_atom_pairs)
            bond_length = round(complete_cohp.bonds[label_list[0]]["length"], 2)
            plot_label = f"{label_list[0]}:  {orb_atom_pairs_str} ({bond_length} \u00c5)"

        return plot_label


class IcohpDistancePlotter:
    """
    Plotter to generate ICOHP or ICOBI or ICOOP vs bond lengths plots.

    :param are_coops: Bool indicating that populations are ICOOPs, not ICOHPs.
                    Defaults to False for COHPs.
    :param are_cobis: Bool indicating that populations are ICOBIs, not ICOHPs.
            Defaults to False for COHPs.
    """

    def __init__(self, are_coops: bool = False, are_cobis: bool = False):
        """Initialize ICOHPs or ICOBI or ICOOP vs bond lengths plotter."""
        self.are_coops = are_coops
        self.are_cobis = are_cobis
        self._icohps = {}  # type: ignore

    def add_icohps(self, label: str, icohpcollection: IcohpCollection):
        """
        Add ICOHPs or ICOBIs or ICOOPS for plotting.

        :param label: Label for the ICOHPs. Must be unique.
        :param icohpcollection: IcohpCollection object.

        """
        icohps = []
        bond_len = []
        atom_pairs = []
        orb_data = {}  # type: ignore
        for indx, bond_label in enumerate(icohpcollection._list_labels):
            orb_data.update({bond_label: {}})
            for k, v in icohpcollection._list_orb_icohp[indx].items():
                orb_data[bond_label].update({k: sum(v["icohp"].values())})
            icohps.append(sum(icohpcollection._list_icohp[indx].values()))
            bond_len.append(icohpcollection._list_length[indx])
            atom1 = icohpcollection._list_atom1[indx]
            atom2 = icohpcollection._list_atom2[indx]
            atom_pairs.append(atom1 + "-" + atom2)

        self._icohps[label] = {
            "atom_pairs": atom_pairs,
            "bond_labels": icohpcollection._list_labels,
            "icohps": icohps,
            "bond_lengths": bond_len,
            "orb_data": orb_data,
        }

    def get_plot(
        self,
        ax: mpl.axes.Axes | None = None,
        marker_size: float = 50,
        marker_style: str = "o",
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        plot_negative: bool = True,
    ):
        """
        Get a matplotlib plot showing the COHP or COBI or COOP with respect to bond lengths.

        :param ax: Existing Matplotlib Axes object to plot to.
        :param marker_size: sets the size of markers in scatter plots
        :param marker_style: sets type of marker used in plot
        :param xlim: Specifies the x-axis limits. Defaults to None for
            automatic determination.
        :param ylim: Specifies the y-axis limits. Defaults to None for
            automatic determination.
        :param plot_negative: Will plot -1*ICOHPs. Works only for ICOHPs

        Returns:
            A matplotlib object.
        """
        if self.are_coops and not self.are_cobis:
            cohp_label = "ICOOP"
        elif self.are_cobis and not self.are_coops:
            cohp_label = "ICOBI"
        elif self.are_cobis and self.are_coops:
            raise ValueError("Plot data should not contain ICOBI and ICOOP data at same time")
        else:
            cohp_label = "ICOHP (eV)"

        if ax is None:
            _, ax = plt.subplots()

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        ax.set_ylabel(cohp_label)
        ax.set_xlabel("Bond lengths (\u00c5)")

        for label, data in self._icohps.items():
            x = data["bond_lengths"]
            if plot_negative and cohp_label == "ICOHP (eV)":
                ax.set_ylabel("$-$" + cohp_label)
                y = [-1 * icohp for icohp in data["icohps"]]
            else:
                y = data["icohps"]

            ax.scatter(x, y, s=marker_size, marker=marker_style, label=label)

        return plt
