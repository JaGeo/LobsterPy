# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
Here classes and functions to plot Lobster outputs are provided
"""

from typing import Any, Tuple, List, Dict

from itertools import cycle
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.electronic_structure.plotter import CohpPlotter
import plotly.graph_objs as go
from lobsterpy.cohp.analyze import Analysis
from lobsterpy.plotting import layout_dicts as ld


base_style = resource_filename("lobsterpy.plotting", "lobsterpy_base.mplstyle")


def get_style_list(
    no_base_style: bool = False,
    styles: "List[str | Dict[str, Any]] | None" = None,
    **kwargs,
) -> "List[str | Dict[str, Any]]":
    """Get *args for matplotlib.style from user input

    Args:
        no_base_style: If true, do not include lobsterpy_base.mplstyle
        styles: User-requested styles. These can be paths to mplstyle files,
                the names of known (matplotlib-supplied) styles,
                or dicts of rcParam options.

    Remaining kwargs are collected as a dict and take highest priority.
    """
    if no_base_style:
        base = []
    else:
        base = [base_style]

    if styles is None:
        styles = []

    return base + styles + [kwargs]


class PlainCohpPlotter(CohpPlotter):
    """
    Modified Pymatgen CohpPlotter with styling removed

    This allows the styling to be manipulated more easily using matplotlib
    style sheets.
    """

    def get_plot(
        self,
        ax: "matplotlib.axes.Axes | None" = None,
        xlim: "Tuple[float, float] | None" = None,
        ylim: "Tuple[float, float] | None" = None,
        plot_negative: "bool | None" = None,
        integrated: bool = False,
        invert_axes: bool = True,
        sigma: "float | None" = None,
    ):
        """
        Get a matplotlib plot showing the COHP.

        Args:
            ax: Existing Matplotlib Axes object to plot to.
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
            plot_negative: It is common to plot -COHP(E) so that the
                sign means the same for COOPs and COHPs. Defaults to None
                for automatic determination: If are_coops is True, this
                will be set to False, else it will be set to True.
            integrated: Switch to plot ICOHPs. Defaults to False.
            invert_axes: Put the energies onto the y-axis, which is
                common in chemistry.
            sigma: Standard deviation of Gaussian broadening applied to
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
            raise ValueError(
                "Plot data should not contain COBI and COOP data at same time"
            )
        else:
            cohp_label = "COHP"

        if plot_negative is None:
            plot_negative = (not self.are_coops) and (not self.are_cobis)

        if integrated:
            cohp_label = "I" + cohp_label + " (eV)"

        if plot_negative:
            cohp_label = "$-$" + cohp_label

        if self.zero_at_efermi:
            energy_label = "$E - E_f$ (eV)"
        else:
            energy_label = "$E$ (eV)"

        colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
        ncolors = len(colors)

        if ax is None:
            _, ax = plt.subplots()

        allpts = []
        keys = self._cohps.keys()
        for i, key in enumerate(keys):
            energies = self._cohps[key]["energies"]
            if not integrated:
                populations = self._cohps[key]["COHP"]
            else:
                populations = self._cohps[key]["ICOHP"]
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
            "color": matplotlib.rcParams["grid.color"],
            "linewidth": matplotlib.rcParams["grid.linewidth"],
            "linestyle": matplotlib.rcParams["grid.linestyle"],
            "alpha": matplotlib.rcParams["grid.alpha"],
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
        """Broaden the spectrum with a given standard deviation

        The population is convolved with a normalised Gaussian kernel. This
        requires the energy grid to be regularly-spaced.

        Args:
            energies: Regularly-spaced energy series
            population: Population data for broadening
            sigma: Standard deviation for Gaussian broadening. If sigma is None
                then the input data is returned without any processing.
            cutoff: Range cutoff for broadening kernel, as a multiple of sigma.

        Return:
            Broadened population
        """
        from scipy.signal import convolve
        from scipy.stats import norm

        if sigma is None:
            return population

        spacing = np.mean(np.diff(energies))
        if not np.allclose(np.diff(energies), spacing, atol=1e-5):
            raise ValueError(
                "Energy grid is not regular, cannot broaden with "
                "discrete convolution."
            )

        # Obtain symmetric mesh for broadening kernel, centered on zero
        kernel_x = np.arange(0, cutoff * sigma + 0.5 * spacing, spacing)
        kernel_x = np.concatenate([-kernel_x[-1:1:-1], kernel_x])

        kernel = norm.pdf(kernel_x, scale=sigma)

        return convolve(population, kernel, mode="same") / kernel.sum()


class InteractiveCohpPlotter(CohpPlotter):
    """
    Interactive COHP plotter to view all relevant / multiple COHPs in one figure.
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

    def add_all_relevant_cohps(
        self, analyse: Analysis, suffix: str = "", label_resolved: bool = True
    ) -> None:
        """
        Adds all relevant COHPs from lobsterpy analyse object.

        Args:
            analyse: Analyse object from lobsterpy.
            suffix: Optional addition to LOBSTER label to avoid key conflicts when plotting multiple
            calcs or just for additional legend information.
            label_resolved: bool indicating to obtain label resolved interactive plots for relevant bonds.
            If false, will return summed cohp curves of unique relevant bonds.
        """
        complete_cohp = analyse.chemenv.completecohp

        # extract bond atom pairs and corresponding cohp bond label
        bonds = [[] for _ in range(len(analyse.seq_infos_bonds))]  # type: ignore
        labels = [[] for _ in range(len(analyse.seq_infos_bonds))]  # type: ignore
        for inx, bond_info in enumerate(analyse.seq_infos_bonds):
            for ixx, val in enumerate(bond_info.atoms):
                label_srt = sorted(val.copy())
                bonds[inx].append(
                    analyse.structure.sites[bond_info.central_isites[0]].species_string
                    + str(bond_info.central_isites[0] + 1)
                    + ": "
                    + label_srt[0].strip("0123456789")
                    + "-"
                    + label_srt[1].strip("0123456789")
                )
                labels[inx].append(bond_info.labels[ixx])

        # create a dict seperating the unique atom pairs for each site and corresponding cohp bond label
        plot_data = {}
        for indx, atom_pairs in enumerate(bonds):
            search_items = set(atom_pairs)
            for item in search_items:
                indices = [i for i, x in enumerate(atom_pairs) if x == item]
                filtered_bond_label_list = [labels[indx][i] for i in indices]
                plot_data.update({item: filtered_bond_label_list})

        if "All" not in self._cohps:
            self._cohps["All"] = {}

        # iterate and extract the data to be plotted from cohp objects
        for bond_key, labels in plot_data.items():
            count = len(labels)
            label_with_count = self._insert_number_of_bonds_in_label(
                label=bond_key, character=":", number_of_bonds=count
            )
            if (
                label_resolved
            ):  # will add cohp data for each relevant bond label iteratively
                self._cohps[label_with_count + suffix] = {}
                for label in labels:
                    cohp = complete_cohp.get_cohp_by_label(label)
                    energies = (
                        cohp.energies - cohp.efermi
                        if self.zero_at_efermi
                        else cohp.energies
                    )
                    outer_key = label_with_count + suffix
                    struct = analyse.structure
                    atom_pairs = []
                    for site in complete_cohp.bonds[label]["sites"]:
                        atom = site.species_string + str(struct.sites.index(site) + 1)
                        atom_pairs.append(atom)
                    key = "{}: {} ({} \u00c5)".format(
                        label,
                        "-".join(atom_pairs),
                        str(round(complete_cohp.bonds[label]["length"], 2)),
                    )
                    self._cohps[outer_key].update(
                        {
                            key: {
                                "energies": energies,
                                "COHP": cohp.get_cohp(),
                                "ICOHP": cohp.get_icohp(),
                                "efermi": cohp.efermi,
                            }
                        }
                    )

                    key = key + suffix
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

            else:
                # add summed cohps for each relevant bond sites
                cohp = complete_cohp.get_summed_cohp_by_label_list(label_list=labels)
                energies = (
                    cohp.energies - cohp.efermi
                    if self.zero_at_efermi
                    else cohp.energies
                )
                key = label_with_count + suffix
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

    def add_cohps_by_lobster_label(
        self, analyse: Analysis, label_list: list, suffix: str = ""
    ):
        """
        Adds COHPs explicitly specified in label list.

        Args:
            analyse: Analyse object from lobsterpy.
            label_list: List of COHP labels as from LOBSTER.
            suffix: Optional addition to LOBSTER label to avoid key
                conflicts when plotting multiple calcs or just for additional legend information.
        """
        complete_cohp = analyse.chemenv.completecohp

        if "All" in self._cohps:
            pass
        else:
            self._cohps["All"] = {}

        for label in label_list:
            atom1 = complete_cohp.bonds[label]["sites"][0].species_string
            atom2 = complete_cohp.bonds[label]["sites"][1].species_string
            sorted_label = sorted([atom1, atom2])
            cohp = complete_cohp.get_cohp_by_label(label)
            energies = (
                cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies
            )
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

    def add_cohps_from_plot_data(self, plot_data_dict: dict, suffix: str = ""):
        """
        Adds all relevant COHPs for specified bond type from lobster lightweight json.gz file

        Args:
            plot_data_dict: Lobsterpy plot data dict
            suffix: Optional addition to LOBSTER label to avoid key
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
                except TypeError:
                    raise ValueError(
                        "The data provided could not be converted to cohp object.Please recheck the input data"
                    )

        if "All" not in self._cohps:
            self._cohps["All"] = {}

        for bond_key, cohps in plot_data_dict.items():
            energies = (
                cohps.energies - cohps.efermi if self.zero_at_efermi else cohps.energies
            )
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
        xlim=None,
        rangeslider=False,
        ylim=None,
        plot_negative=None,
        integrated=False,
        invert_axes=True,
        sigma=None,
        colors=None,
    ):
        """
        Get an interactive plotly figure showing the COHPs.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            rangeslider: Adds a plotly.graph_objs.layout.xaxis.Rangeslider
                object to figure to allow easy manipulation of x-axis..
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
            plot_negative: It is common to plot -COHP(E) so that the
                sign means the same for COOPs and COHPs. Defaults to None
                for automatic determination: If are_coops is True, this
                will be set to False, else it will be set to True.
            integrated: Switch to plot ICOHPs. Defaults to False.
            invert_axes: Put the energies onto the y-axis, which is
                common in chemistry.
            sigma: Standard deviation of Gaussian broadening applied to
                population data. If this is unset (None) no broadening will be
                added.
            colors: list of hex color codes to be used in plot

        Returns:
            A  plotly.graph_objects.Figure object.
        """
        if self.are_coops and not self.are_cobis:
            cohp_label = "COOP"
        elif self.are_cobis and not self.are_coops:
            cohp_label = "COBI"
        elif self.are_cobis and self.are_coops:
            raise ValueError(
                "Plot data should not contain COBI and COOP data at same time"
            )
        else:
            cohp_label = "COHP"

        if plot_negative is None:
            plot_negative = (not self.are_coops) and (not self.are_cobis)

        if integrated:
            cohp_label = "I" + cohp_label + " (eV)"

        if plot_negative:
            cohp_label = "\u2212" + cohp_label

        if self.zero_at_efermi:
            energy_label = "$E - E_f$ (eV)"
        else:
            energy_label = "$E$ (eV)"

        # Setting up repeating color scheme (same as for matplotlib plots in .mplstyle)
        if colors is None:
            palette = InteractiveCohpPlotter.COLOR_PALETTE
        else:
            palette = colors

        pal_iter = cycle(palette)

        traces = {}
        for k, v in self._cohps.items():
            traces.update({k: []})
            for label, val in v.items():
                population_key = v[label]["ICOHP"] if integrated else v[label]["COHP"]
                band_color = next(pal_iter)
                for spin in [Spin.up, Spin.down]:
                    if spin in population_key:
                        population = (
                            [-i for i in population_key[spin]]
                            if plot_negative
                            else population_key[spin]
                        )
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
            else go.layout.XAxis(
                title=energy_label, rangeslider={"visible": rangeslider}
            )
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
        for _, val_trace in traces.items():
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
                            for selected_group in traces.keys()
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
    def _insert_number_of_bonds_in_label(
        label: str, character: str, number_of_bonds: int
    ) -> str:
        """
        Adds number of bonds to bond label.
        For example : for input label 'Ba1: Ba-Ti', character ':', number_of_bonds: 3,
        Will return 'Ba1: 3 x Ba-Ti'

        Args:
            label: bond label to which number of bonds needs to be inserted
            character: string character where number of bonds needs to be inserted
            number_of_bonds: number of bonds corresponding to the label

        Returns:
             bond label with number of bonds inserted
        """
        return label.replace(character, f"{character} {number_of_bonds} x", 1)
