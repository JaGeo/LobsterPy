# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
Here classes and functions to plot Lobster outputs are provided
"""

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.electronic_structure.core import Spin

base_style = resource_filename("lobsterpy.plotting", "lobsterpy_base.mplstyle")


def get_style_list(
    no_base_style: bool = False,
    styles: Optional[List[Union[str, Dict[str, Any]]]] = None,
    **kwargs
) -> List[Union[str, Dict[str, Any]]]:
    """Get *args for matplotlib.style from user input

    Args:
        no_base_style: If true, do not include lobsterpy_base.mplstyle
        styles: User-requested styles. These can be paths to mplstyle files,
                the names of known (matplotlib-supplied) styles,
                or dicts of rcParam options.

    Remaining kwargs are collected as a dict and take highest priority.

    """

    if no_base_style:
        base = []  # type: List[Union[str, Dict[str, Any]]]
    else:
        base = [base_style]

    if styles is None:
        styles = []

    return base + styles + [kwargs]


class PlainCohpPlotter(CohpPlotter):
    """Modified Pymatgen CohpPlotter with styling removed

    This allows the styling to be manipulated more easily using matplotlib
    style sheets."""

    def get_plot(
        self,
        ax: Optional[matplotlib.axes.Axes] = None,
        xlim: Optional[Tuple[float, float]] = None,
        ylim: Optional[Tuple[float, float]] = None,
        plot_negative: Optional[bool] = None,
        integrated: bool = False,
        invert_axes: bool = True,
        sigma: Optional[float] = None,
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
        if self.are_coops:
            cohp_label = "COOP"
        elif self.are_cobis:
            cohp_label = "COBI"
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
