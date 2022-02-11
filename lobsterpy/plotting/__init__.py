from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from pymatgen.electronic_structure.plotter import CohpPlotter

base_style = resource_filename('lobsterpy.plotting', 'lobsterpy_base.mplstyle')


class PlainCohpPlotter(CohpPlotter):
    """Modified Pymatgen CohpPlotter with styling removed

    This allows the styling to be manipulated more easily using matplotlib
    style sheets."""

    def get_plot(
        self,
        ax=None,
        xlim=None,
        ylim=None,
        plot_negative=None,
        integrated=False,
        invert_axes=True,
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
            cohp_label = "-" + cohp_label

        if self.zero_at_efermi:
            energy_label = "$E - E_f$ (eV)"
        else:
            energy_label = "$E$ (eV)"

        colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
        ncolors = len(colors)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

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
                    else:
                        x = energies
                        y = -populations[spin] if plot_negative else populations[spin]
                    allpts.extend(list(zip(x, y)))
                    if spin == Spin.up:
                        plt.plot(
                            x,
                            y,
                            color=colors[i % ncolors],
                            linestyle="-",
                            label=str(key),
                            linewidth=3,
                        )
                    else:
                        plt.plot(x, y, color=colors[i % ncolors], linestyle="--", linewidth=3)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        xlim = plt.xlim()
        ylim = plt.ylim()
        if not invert_axes:
            plt.plot(xlim, [0, 0], "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot([0, 0], ylim, "k--", linewidth=2)
            else:
                plt.plot(
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    ylim,
                    color=colors[i % ncolors],
                    linestyle="--",
                    linewidth=2,
                )
        else:
            plt.plot([0, 0], ylim, "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot(xlim, [0, 0], "k--", linewidth=2)
            else:
                plt.plot(
                    xlim,
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    color=colors[i % ncolors],
                    linestyle="--",
                    linewidth=2,
                )

        if invert_axes:
            plt.xlabel(cohp_label)
            plt.ylabel(energy_label)
        else:
            plt.xlabel(energy_label)
            plt.ylabel(cohp_label)

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt
