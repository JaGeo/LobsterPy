from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from plotly.io import read_json
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.electronic_structure.core import Spin

from lobsterpy.cohp.describe import Description
from lobsterpy.plotting import (
    IcohpDistancePlotter,
    InteractiveCohpPlotter,
    PlainCohpPlotter,
    PlainDosPlotter,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestInteractiveCohpPlotter:
    def test_add_all_relevant_cohps_nacl(self, plot_analyse_nacl):
        iplotter = InteractiveCohpPlotter(zero_at_efermi=False)

        iplotter.add_all_relevant_cohps(analyse=plot_analyse_nacl, label_resolved=False, suffix="")
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 1

        fig = iplotter.get_plot(invert_axes=False)
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_NaCl.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(list(og_trace.y))
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_all_relevant_cohps_nacl_cobi_orb(self, plot_analyse_nacl_cobi_orb):
        iplotter = InteractiveCohpPlotter(are_cobis=True)

        iplotter.add_all_relevant_cohps(
            analyse=plot_analyse_nacl_cobi_orb,
            label_resolved=False,
            suffix="",
            orbital_resolved=True,
        )
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 3

        fig = iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_NaCl_cobi_orb_res.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_all_relevant_cohps_batio3_orb(self, plot_analyse_batio3_orb):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_all_relevant_cohps(
            analyse=plot_analyse_batio3_orb,
            label_resolved=True,
            orbital_resolved=True,
            suffix="",
        )
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 6

        fig = iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_BaTiO3_orb_label_resol.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_all_relevant_cohps_cdf(self, plot_analyse_cdf_orb):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_all_relevant_cohps(
            analyse=plot_analyse_cdf_orb,
            label_resolved=False,
            orbital_resolved=True,
            suffix="",
        )
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 3

        fig = iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_CdF_orb_resol.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_all_relevant_cohps_nacl_cobi(self, plot_analyse_nacl_cobi):
        iplotter = InteractiveCohpPlotter(zero_at_efermi=False, are_cobis=True)

        iplotter.add_all_relevant_cohps(analyse=plot_analyse_nacl_cobi, label_resolved=False, suffix="")
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 1

        fig = iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_NaCl_cobi.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_all_relevant_cohps_k3sb(self, plot_analyse_k3sb):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_all_relevant_cohps(analyse=plot_analyse_k3sb, label_resolved=True, suffix="")
        assert "All" in iplotter._cohps
        assert "K1: 8 x K-K" in iplotter._cohps
        assert "K1: 6 x K-Sb" in iplotter._cohps
        assert "K2: 4 x K-Sb" in iplotter._cohps
        assert "K2: 10 x K-K" in iplotter._cohps
        assert "K2: 10 x K-K" in iplotter._cohps
        assert "Sb4: 14 x K-Sb" in iplotter._cohps
        assert len(iplotter._cohps) == 6

        fig = iplotter.get_plot(sigma=0.3, xlim=[-5, 5], ylim=[-10, 10])
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_K3Sb.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout.xaxis == ref_fig.layout.xaxis
        assert fig.layout.yaxis == ref_fig.layout.yaxis
        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_cohps_by_lobster_label_nacl(self, plot_analyse_nacl):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_cohps_by_lobster_label(analyse=plot_analyse_nacl, label_list=["5", "10", "15"], suffix="")
        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 1

        fig = iplotter.get_plot(integrated=True)
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_NaCl_label.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout
        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_cohps_from_plot_data(self, plot_analyse_nasi):
        des = Description(analysis_object=plot_analyse_nasi)

        fig = des.plot_interactive_cohps(hide=True)
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/analyse_NaSi.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_cohps_from_plot_data_json(self, lobsterpy_plot_data):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_cohps_from_plot_data(plot_data_dict=lobsterpy_plot_data, suffix="")

        assert "All" in iplotter._cohps
        assert len(iplotter._cohps) == 1

        fig = iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "test_data/interactive_plotter_ref/fig_mp8818.json",
            engine="json",
        )
        assert len(fig.data) == len(ref_fig.data)
        assert fig.layout == ref_fig.layout

        og_name, og_line, og_visible, og_x, og_y = [], [], [], [], []
        ref_name, ref_line, ref_visible, ref_x, ref_y = [], [], [], [], []
        for og_trace, ref_trace in zip(fig.data, ref_fig.data):
            og_name.append(og_trace.name)
            ref_name.append(ref_trace.name)
            og_line.append(og_trace.line)
            ref_line.append(ref_trace.line)
            og_visible.append(og_trace.visible)
            ref_visible.append(ref_trace.visible)
            og_x.append(list(og_trace.x))
            ref_x.append(list(ref_trace.x))
            og_y.append(og_trace.y.tolist())
            ref_y.append(list(ref_trace.y))

        assert sorted(og_name) == sorted(ref_name)
        assert og_line == ref_line
        assert og_visible == ref_visible

        # use numpy testing for array equality
        np.testing.assert_array_almost_equal(sorted(og_x), sorted(ref_x), decimal=4)
        np.testing.assert_array_almost_equal(sorted(og_y), sorted(ref_y), decimal=4)

    def test_add_cohp_nacl(self, plot_analyse_nacl):
        iplotter = InteractiveCohpPlotter()

        for _, (ication, labels, cohps) in enumerate(
            zip(
                plot_analyse_nacl.seq_ineq_ions,
                plot_analyse_nacl.seq_labels_cohps,
                plot_analyse_nacl.seq_cohps,
            )
        ):
            namecation = str(plot_analyse_nacl.structure[ication].specie)
            for label, cohp in zip(labels, cohps):
                if label is not None:
                    iplotter.add_cohp(namecation + str(ication + 1) + ": " + label, cohp)

        fig = iplotter.get_plot()

        assert len(fig.data) == 1

    def test_plot_colors(self, plot_analyse_k3sb):
        iplotter = InteractiveCohpPlotter()

        iplotter.add_all_relevant_cohps(analyse=plot_analyse_k3sb, label_resolved=True, suffix="")

        Fig_ref = iplotter.get_plot()

        color = ["#00FFFF", "#008080", "#00008B", "#808000"]
        # custom color plot
        Fig_cust_col = iplotter.get_plot(colors=color)

        for ref, cust in zip(Fig_ref.data, Fig_cust_col.data):
            assert ref.line.color != cust.line.color

    def test_plot_labels(self):
        # plain cohp plotter
        plotter = PlainCohpPlotter(are_cobis=True)
        fig = plotter.get_plot().gca()

        assert fig.get_xlabel() == "COBI"

        plotter = PlainCohpPlotter(are_coops=True)
        fig = plotter.get_plot().gca()

        assert fig.get_xlabel() == "COOP"

        plotter = PlainCohpPlotter()
        fig = plotter.get_plot().gca()

        assert fig.get_xlabel() == "$-$COHP (eV)"

        # interactive plotter
        iplotter = InteractiveCohpPlotter(are_cobis=True)
        fig = iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "COBI"

        iplotter = InteractiveCohpPlotter(are_coops=True)
        fig = iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "COOP"

        iplotter = InteractiveCohpPlotter()
        fig = iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "−COHP (eV)"  # noqa: RUF001
        assert fig.layout.yaxis["title"]["text"] == "$E - E_f \\text{ (eV)}$"

    def test_plaincohp_plotter_options(self, lobsterpy_plot_data):
        plotter = PlainCohpPlotter(zero_at_efermi=False)

        for label, cohp in lobsterpy_plot_data.items():
            cohp_obj = Cohp.from_dict(cohp)
            plotter.add_cohp(label=label, cohp=cohp_obj)

        fig = plotter.get_plot(integrated=True, xlim=(-5, 2), ylim=(-4, 4), invert_axes=False).gca()
        assert fig.get_ylabel() == "$-$ICOHP (eV)"
        assert fig.get_xlabel() == "$E$ (eV)"


class TestIcohpDistancePlotter:
    def test_icohp_plotter_labels(self, icohplist_nacl, icooplist_nacl, icobilist_nacl):
        icohp_plotter = IcohpDistancePlotter()
        icohp_plotter.add_icohps(label="NaCl", icohpcollection=icohplist_nacl.icohpcollection)
        fig = icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "$-$" + "ICOHP (eV)"

        icohp_plotter = IcohpDistancePlotter()
        icohp_plotter.add_icohps(label="NaCl_icohp", icohpcollection=icohplist_nacl.icohpcollection)
        fig = icohp_plotter.get_plot(plot_negative=True).gca()
        assert fig.get_ylabel() == "$-$" + "ICOHP (eV)"

        icohp_plotter = IcohpDistancePlotter(are_cobis=True)
        icohp_plotter.add_icohps(label="NaCl_icobi", icohpcollection=icobilist_nacl.icohpcollection)
        fig = icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "ICOBI"

        icohp_plotter = IcohpDistancePlotter(are_coops=True)
        icohp_plotter.add_icohps(label="NaCl_icoop", icohpcollection=icooplist_nacl.icohpcollection)
        fig = icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "ICOOP"
        assert fig.get_xlabel() == "Bond lengths (Å)"

    def test_plot_data(self, icohplist_nacl, icobilist_nacl, icooplist_nacl):
        icohp_plotter = IcohpDistancePlotter()
        icohp_plotter.add_icohps(label="NaCl", icohpcollection=icohplist_nacl.icohpcollection)

        ref_xdata = icohplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in icohplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(abs(sum(ydata.values())))  # get absolute icohp values as in plots

        fig_xydata = icohp_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = icohp_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert fig_x_lims == [0, 4]
        assert fig_y_lims == [0, 6]

        # icobi

        icobi_plotter = IcohpDistancePlotter(are_cobis=True)
        icobi_plotter.add_icohps(label="NaCl", icohpcollection=icobilist_nacl.icohpcollection)

        ref_xdata = icobilist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in icobilist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(sum(ydata.values()))  # get absolute icohp values as in plots

        fig_xydata = icobi_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = icobi_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert fig_x_lims == [0, 4]
        assert fig_y_lims == [0, 6]

        # icoop

        icoop_plotter = IcohpDistancePlotter(are_coops=True)
        icoop_plotter.add_icohps(label="NaCl", icohpcollection=icooplist_nacl.icohpcollection)

        ref_xdata = icooplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in icooplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(sum(ydata.values()))  # get absolute icohp values as in plots

        fig_xydata = icoop_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = icoop_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert fig_x_lims == [0, 4]
        assert fig_y_lims == [0, 6]

        # test for colors in icoxx plotter
        ax_icoop_c = icoop_plotter.get_plot(xlim=(0, 4), ylim=(0, 6), color_interactions=True).gca()
        ax_icohp_c = icohp_plotter.get_plot(xlim=(0, 4), ylim=(0, 6), color_interactions=True).gca()
        ax_icobi_c = icobi_plotter.get_plot(xlim=(0, 4), ylim=(0, 6), color_interactions=True).gca()
        handles_icoop, labels_icoop = ax_icoop_c.get_legend_handles_labels()
        handles_icohp, labels_icohp = ax_icohp_c.get_legend_handles_labels()
        handles_icobi, labels_icobi = ax_icobi_c.get_legend_handles_labels()
        assert len(set(labels_icoop)) == 3
        assert len(set(labels_icobi)) == 3
        assert len(set(labels_icohp)) == 3


class TestPlotterExceptions:
    def test_plotter_exception(self, plot_analyse_nasi):
        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            iplotter = InteractiveCohpPlotter()

            data = {"N4: 1 x N-N": []}

            iplotter.add_cohps_from_plot_data(plot_data_dict=data, suffix="")

        assert str(err.value) == "The data provided could not be converted to cohp object.Please recheck the input data"

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            iplotter = InteractiveCohpPlotter(are_cobis=True, are_coops=True)

            _ = iplotter.get_plot()

        assert str(err.value) == "Plot data should not contain COBI and COOP data at same time"

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            plotter = PlainCohpPlotter(are_cobis=True, are_coops=True)

            _ = plotter.get_plot()

        assert str(err.value) == "Plot data should not contain COBI and COOP data at same time"

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            iplotter = InteractiveCohpPlotter()
            iplotter.add_cohp(
                label="10",
                cohp=plot_analyse_nasi.chemenv.completecohp.get_cohp_by_label("10"),
            )
            plotter.add_cohp(
                label="10",
                cohp=plot_analyse_nasi.chemenv.completecohp.get_cohp_by_label("20"),
            )

            assert (
                str(err.value) == "Please use another label to add the COHP, provided label already exists "
                "in the plot data, which will result in overwriting the existing COHP data."
            )

        with pytest.raises(Exception) as err:  # noqa: PT012, PT011
            icohp_plotter = IcohpDistancePlotter(are_coops=True, are_cobis=True)

            _ = icohp_plotter.get_plot()

            assert str(err.value) == "Plot data should not contain ICOBI and ICOOP data at same time"


class TestPlainDosPlotter:
    def test_nacl_dos(self, nacl_dos):
        complete_dos_obj = nacl_dos.completedos
        # add and test total non normalized dos data and axis labels in the plot
        dp = PlainDosPlotter(summed=False, stack=False, sigma=None)
        dp.add_dos(dos=complete_dos_obj, label="Total")
        plt = dp.get_plot(invert_axes=False, beta_dashed=True).gcf()

        for energies in plt.axes[0].get_lines()[:-2]:
            plot_en = energies.get_data()[0].tolist()
            ref_en = complete_dos_obj.energies.tolist()
            assert plot_en == ref_en

        for plot_dos, ref_dos in zip(plt.axes[0].get_lines()[:-2], complete_dos_obj.densities.values()):
            dos_plot = [abs(dos) for dos in plot_dos.get_data()[1].tolist()]
            dos_ref = [abs(dos) for dos in ref_dos.tolist()]
            assert dos_plot == dos_ref

        plt_axes = dp.get_plot(invert_axes=False, beta_dashed=True).gca()

        assert plt_axes.get_xlabel() == "Energies (eV)"
        assert plt_axes.get_ylabel() == "Density of states (states/eV)"

        # add and test total normalized dos data and axis labels in the plot

        complete_dos_obj_norm = nacl_dos.completedos.get_normalized()
        dp = PlainDosPlotter(summed=False, stack=False, sigma=None)
        dp.add_dos(dos=complete_dos_obj_norm, label="Total")
        plt = dp.get_plot(invert_axes=False, beta_dashed=True).gcf()

        for energies in plt.axes[0].get_lines()[:-2]:
            plot_en = energies.get_data()[0].tolist()
            ref_en = complete_dos_obj_norm.energies.tolist()
            assert plot_en == ref_en

        for plot_dos, ref_dos in zip(plt.axes[0].get_lines()[:-2], complete_dos_obj_norm.densities.values()):
            dos_plot = [abs(dos) for dos in plot_dos.get_data()[1].tolist()]
            dos_ref = [abs(dos) for dos in ref_dos.tolist()]
            assert dos_plot == dos_ref

        plt_axes = dp.get_plot(invert_axes=False, beta_dashed=True).gca()

        assert plt_axes.get_xlabel() == "Energies (eV)"
        assert plt_axes.get_ylabel() == "Density of states (states/eV/Å³)"

    def test_k3sb_dos(self, k3sb_dos):
        complete_dos_obj = k3sb_dos.completedos
        # add and test total non normalized dos data and axis labels in the plot
        dp = PlainDosPlotter(summed=True, stack=False, sigma=None)
        dp.add_dos(dos=complete_dos_obj, label="Total")
        plt = dp.get_plot(invert_axes=True, beta_dashed=True).gcf()

        for energies in plt.axes[0].get_lines()[:1]:
            plot_en = energies.get_data()[1].tolist()
            ref_en = complete_dos_obj.energies.tolist()
            assert plot_en == ref_en

        for plot_dos in plt.axes[0].get_lines()[:1]:
            dos_plot = [abs(dos) for dos in plot_dos.get_data()[0].tolist()]
            dos_ref = [abs(dos) for dos in complete_dos_obj.get_densities().tolist()]
            assert dos_plot == dos_ref

        plt_axes = dp.get_plot(invert_axes=True, beta_dashed=True).gca()

        assert plt_axes.get_ylabel() == "Energies (eV)"
        assert plt_axes.get_xlabel() == "Density of states (states/eV)"

        # add and test total non-normalized smeared dos data and axis labels in the plot
        dp = PlainDosPlotter(summed=True, stack=False, sigma=0.1)
        dp.add_dos(dos=complete_dos_obj, label="Total")
        plt = dp.get_plot(invert_axes=False, beta_dashed=True).gcf()

        for energies in plt.axes[0].get_lines()[:1]:
            plot_en = energies.get_data()[0].tolist()
            ref_en = complete_dos_obj.energies.tolist()
            assert plot_en == ref_en

        for plot_dos in plt.axes[0].get_lines()[:1]:
            dos_plot = [abs(dos) for dos in plot_dos.get_data()[1].tolist()]
            dos_ref = [abs(dos) for dos in sum(complete_dos_obj.get_smeared_densities(sigma=0.1).values()).tolist()]
            assert dos_plot == dos_ref

        plt_axes = dp.get_plot(invert_axes=False, beta_dashed=True).gca()

        assert plt_axes.get_xlabel() == "Energies (eV)"
        assert plt_axes.get_ylabel() == "Density of states (states/eV)"

    def test_add_site_orbital_dos(self, nacl_dos, k3sb_dos):
        # test for all orbitals plot at the site
        complete_dos_k3sb = k3sb_dos.completedos
        dp = PlainDosPlotter(summed=True, stack=False, sigma=0.01)
        dp.add_site_orbital_dos(dos=complete_dos_k3sb, site_index=0, orbital="all")

        plt = dp.get_plot().gcf()

        site = complete_dos_k3sb.structure.sites[0]
        avail_orbs = list(complete_dos_k3sb.pdos[site])

        for data in plt.axes[0].get_lines()[:-2]:
            if data.get_label().split(":")[-1].strip() in avail_orbs:
                dos_obj = complete_dos_k3sb.get_site_orbital_dos(
                    site=site, orbital=data.get_label().split(":")[-1].strip()
                )
                plot_en = data.get_data()[0].tolist()
                plot_dos = data.get_data()[1].tolist()
                ref_en = dos_obj.energies.tolist()
                ref_dos = dos_obj.get_smeared_densities(0.01)
                ref_dos_add = ref_dos[Spin.up] + ref_dos[Spin.down]
                np.testing.assert_array_almost_equal(sorted(plot_en), sorted(ref_en), decimal=4)
                np.testing.assert_array_almost_equal(sorted(plot_dos), sorted(ref_dos_add), decimal=4)
            else:
                raise Exception("Plot data does not match expected output")

        # test for specific orbital
        complete_dos_nacl = nacl_dos.completedos
        dp = PlainDosPlotter(summed=True, stack=False, sigma=0.01)
        dp.add_site_orbital_dos(dos=complete_dos_nacl, site_index=0, orbital="3s")

        plt = dp.get_plot().gcf()

        site2 = complete_dos_nacl.structure.sites[0]

        for data in plt.axes[0].get_lines()[:-2]:
            dos_obj = complete_dos_nacl.get_site_orbital_dos(site=site2, orbital="3s")
            plot_en = data.get_data()[0].tolist()
            plot_dos = data.get_data()[1].tolist()
            ref_en = dos_obj.energies.tolist()
            ref_dos = dos_obj.get_smeared_densities(0.01)
            ref_dos_add = ref_dos[Spin.up] + ref_dos[Spin.down]
            np.testing.assert_array_almost_equal(sorted(plot_en), sorted(ref_en), decimal=4)
            np.testing.assert_array_almost_equal(sorted(plot_dos), sorted(ref_dos_add), decimal=4)

    def test_dos_plotter_exceptions(self, nacl_dos):
        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            dp = PlainDosPlotter(summed=True, stack=False, sigma=None)

            _ = dp.add_site_orbital_dos(site_index=0, orbital="5_s", dos=nacl_dos.completedos)

        assert (
            str(err.value)
            == "Requested orbital is not available for this site, available orbitals are 3s, 2p_y, 2p_z, 2p_x"
        )
