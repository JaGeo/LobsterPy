from __future__ import annotations

import unittest
import gzip
import json
from pathlib import Path
from plotly.io import read_json
from pymatgen.io.lobster import Icohplist
from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from lobsterpy.plotting import (
    PlainCohpPlotter,
    InteractiveCohpPlotter,
    IcohpDistancePlotter,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


class InteractiveCohpPlotterTest(unittest.TestCase):
    def setUp(self):
        self.analyse_NaCl = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            summed_spins=False,
        )

        self.analyse_NaSi = Analysis(
            path_to_poscar=TestDir / "TestData/NaSi/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaSi/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaSi/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaSi/CHARGE.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
            summed_spins=True,
        )

        self.analyse_K3Sb = Analysis(
            path_to_poscar=TestDir / "TestData/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/K3Sb/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/K3Sb/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/K3Sb/CHARGE.lobster.gz",
            whichbonds="all",
            cutoff_icohp=0.1,
            summed_spins=False,
        )

        plot_data_file_name = (
            TestDir / "TestData/interactive_plotter_ref/mp-8818.json.gz"
        )

        with gzip.open(plot_data_file_name, "rb") as f:
            data = json.loads(f.read().decode("utf-8"))

        lobsterpy_plot_data = {}
        for item in data:
            lobsterpy_plot_data.update(item)

        self.lobsterpy_plot_data = lobsterpy_plot_data["all_bonds"]["lobsterpy_data"][
            "cohp_plot_data"
        ]

    def test_add_all_relevant_cohps_NaCl(self):
        self.iplotter = InteractiveCohpPlotter(zero_at_efermi=False)

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_NaCl, label_resolved=False, suffix=""
        )
        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 1)

        fig = self.iplotter.get_plot(invert_axes=False)
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_NaCl.json",
            engine="json",
        )
        self.assertEqual(len(fig.data), len(ref_fig.data))
        self.assertEqual(fig.layout, ref_fig.layout)
        for og_trace in fig.data:
            if og_trace in ref_fig.data:
                ref_trace = ref_fig.data[ref_fig.data.index(og_trace)]
                for og_x, og_y, ref_x, ref_y in zip(
                    og_trace.x, og_trace.y, ref_trace.x, ref_trace.y
                ):
                    self.assertAlmostEqual(ref_x, og_x, delta=0.0001)
                    self.assertAlmostEqual(ref_y, og_y, delta=0.0001)
                self.assertEqual(og_trace.name, ref_trace.name)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.visible, ref_trace.visible)

    def test_add_all_relevant_cohps_K3Sb(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_K3Sb, label_resolved=True, suffix=""
        )
        self.assertIn("All", self.iplotter._cohps)
        self.assertIn("K1: 8 x K-K", self.iplotter._cohps)
        self.assertIn("K1: 6 x K-Sb", self.iplotter._cohps)
        self.assertIn("K2: 4 x K-Sb", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("Sb4: 14 x K-Sb", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 6)

        fig = self.iplotter.get_plot(sigma=0.3, xlim=[-5, 5], ylim=[-10, 10])
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_K3Sb.json",
            engine="json",
        )
        self.assertEqual(len(fig.data), len(ref_fig.data))
        self.assertEqual(fig.layout.xaxis, ref_fig.layout.xaxis)
        self.assertEqual(fig.layout.yaxis, ref_fig.layout.yaxis)
        for og_trace in fig.data:
            if og_trace in ref_fig.data:
                ref_trace = ref_fig.data[ref_fig.data.index(og_trace)]
                for og_x, og_y, ref_x, ref_y in zip(
                    og_trace.x, og_trace.y, ref_trace.x, ref_trace.y
                ):
                    self.assertAlmostEqual(ref_x, og_x, delta=0.0001)
                    self.assertAlmostEqual(ref_y, og_y, delta=0.0001)
                self.assertEqual(og_trace.name, ref_trace.name)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.visible, ref_trace.visible)

    def test_add_cohps_by_lobster_label_NaCl(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_cohps_by_lobster_label(
            analyse=self.analyse_NaCl, label_list=["5", "10", "15"], suffix=""
        )
        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 1)

        fig = self.iplotter.get_plot(integrated=True)
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_NaCl_label.json",
            engine="json",
        )
        self.assertEqual(len(fig.data), len(ref_fig.data))
        self.assertEqual(fig.layout, ref_fig.layout)
        for og_trace in fig.data:
            if og_trace in ref_fig.data:
                ref_trace = ref_fig.data[ref_fig.data.index(og_trace)]
                for og_x, og_y, ref_x, ref_y in zip(
                    og_trace.x, og_trace.y, ref_trace.x, ref_trace.y
                ):
                    self.assertAlmostEqual(ref_x, og_x, delta=0.0001)
                    self.assertAlmostEqual(ref_y, og_y, delta=0.0001)
                self.assertEqual(og_trace.name, ref_trace.name)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.visible, ref_trace.visible)

    def test_add_cohps_from_plot_data(self):
        self.des = Description(analysis_object=self.analyse_NaSi)

        fig = self.des.plot_interactive_cohps(hide=True)
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_NaSi.json",
            engine="json",
        )
        self.assertEqual(len(fig.data), len(ref_fig.data))
        self.assertEqual(fig.layout, ref_fig.layout)

        for og_trace in fig.data:
            if og_trace in ref_fig.data:
                ref_trace = ref_fig.data[ref_fig.data.index(og_trace)]
                for og_x, og_y, ref_x, ref_y in zip(
                    og_trace.x, og_trace.y, ref_trace.x, ref_trace.y
                ):
                    self.assertAlmostEqual(ref_x, og_x, delta=0.0001)
                    self.assertAlmostEqual(ref_y, og_y, delta=0.0001)
                self.assertEqual(og_trace.name, ref_trace.name)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.visible, ref_trace.visible)

    def test_add_cohps_from_plot_data_json(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_cohps_from_plot_data(
            plot_data_dict=self.lobsterpy_plot_data, suffix=""
        )

        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 1)

        fig = self.iplotter.get_plot()
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/fig_mp8818.json",
            engine="json",
        )
        self.assertEqual(len(fig.data), len(ref_fig.data))
        self.assertEqual(fig.layout, ref_fig.layout)

        for og_trace in fig.data:
            if og_trace in ref_fig.data:
                ref_trace = ref_fig.data[ref_fig.data.index(og_trace)]
                for og_x, og_y, ref_x, ref_y in zip(
                    og_trace.x, og_trace.y, ref_trace.x, ref_trace.y
                ):
                    self.assertAlmostEqual(ref_x, og_x, delta=0.0001)
                    self.assertAlmostEqual(ref_y, og_y, delta=0.0001)
                self.assertEqual(og_trace.name, ref_trace.name)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.line, ref_trace.line)
                self.assertEqual(og_trace.visible, ref_trace.visible)

    def test_plot_colors(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_K3Sb, label_resolved=True, suffix=""
        )

        Fig_ref = self.iplotter.get_plot()

        color = ["#00FFFF", "#008080", "#00008B", "#808000"]
        # custom color plot
        Fig_cust_col = self.iplotter.get_plot(colors=color)

        for ref, cust in zip(Fig_ref.data, Fig_cust_col.data):
            self.assertNotEqual(ref.line.color, cust.line.color)

    def test_plot_labels(self):
        # plain cohp plotter
        self.plotter = PlainCohpPlotter(are_cobis=True)
        fig = self.plotter.get_plot().gca()

        self.assertEqual(fig.get_xlabel(), "COBI")

        self.plotter = PlainCohpPlotter(are_coops=True)
        fig = self.plotter.get_plot().gca()

        self.assertEqual(fig.get_xlabel(), "COOP")

        self.plotter = PlainCohpPlotter()
        fig = self.plotter.get_plot().gca()

        self.assertEqual(fig.get_xlabel(), "$-$COHP (eV)")

        # interactive plotter
        self.iplotter = InteractiveCohpPlotter(are_cobis=True)
        fig = self.iplotter.get_plot()

        self.assertEqual(fig.layout.xaxis["title"]["text"], "COBI")

        self.iplotter = InteractiveCohpPlotter(are_coops=True)
        fig = self.iplotter.get_plot()

        self.assertEqual(fig.layout.xaxis["title"]["text"], "COOP")

        self.iplotter = InteractiveCohpPlotter()
        fig = self.iplotter.get_plot()

        self.assertEqual(fig.layout.xaxis["title"]["text"], "−COHP (eV)")
        self.assertEqual(fig.layout.yaxis["title"]["text"], "$E - E_f \\text{ (eV)}$")


class IcohpDistancePlotterTest(unittest.TestCase):
    def setUp(self):
        self.icohplist_nacl = Icohplist(
            filename=TestDir / "TestData/NaCl_comp_range/ICOHPLIST.lobster.gz"
        )
        self.icobilist_nacl = Icohplist(
            filename=TestDir / "TestData/NaCl_comp_range/ICOBILIST.lobster.gz",
            are_cobis=True,
        )
        self.icooplist_nacl = Icohplist(
            filename=TestDir / "TestData/NaCl_comp_range/ICOOPLIST.lobster.gz",
            are_coops=True,
        )

    def test_icohp_plotter_labels(self):
        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(
            label="NaCl", icohpcollection=self.icohplist_nacl.icohpcollection
        )
        fig = self.icohp_plotter.get_plot().gca()
        self.assertEqual(fig.get_ylabel(), "$-$" + "ICOHP (eV)")

        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(
            label="NaCl_icohp", icohpcollection=self.icohplist_nacl.icohpcollection
        )
        fig = self.icohp_plotter.get_plot(plot_negative=True).gca()
        self.assertEqual(fig.get_ylabel(), "$-$" + "ICOHP (eV)")

        self.icohp_plotter = IcohpDistancePlotter(are_cobis=True)
        self.icohp_plotter.add_icohps(
            label="NaCl_icobi", icohpcollection=self.icobilist_nacl.icohpcollection
        )
        fig = self.icohp_plotter.get_plot().gca()
        self.assertEqual(fig.get_ylabel(), "ICOBI")

        self.icohp_plotter = IcohpDistancePlotter(are_coops=True)
        self.icohp_plotter.add_icohps(
            label="NaCl_icoop", icohpcollection=self.icooplist_nacl.icohpcollection
        )
        fig = self.icohp_plotter.get_plot().gca()
        self.assertEqual(fig.get_ylabel(), "ICOOP")
        self.assertEqual(fig.get_xlabel(), "Bond lengths (Å)")

    def test_plot_data(self):
        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(
            label="NaCl", icohpcollection=self.icohplist_nacl.icohpcollection
        )

        ref_xdata = self.icohplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icohplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(
                abs(sum(ydata.values()))
            )  # get absolute icohp values as in plots

        fig_xydata = (
            self.icohp_plotter.get_plot()
            .gcf()
            .axes[0]
            .get_children()[0]
            .get_offsets()
            .data
        )

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        self.assertListEqual(ref_xdata, fig_xdata)
        self.assertListEqual(ref_ydata, fig_ydata)

        fig_xydata = self.icohp_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        self.assertListEqual([0, 4], fig_x_lims)
        self.assertListEqual([0, 6], fig_y_lims)

        # icobi

        self.icobi_plotter = IcohpDistancePlotter(are_cobis=True)
        self.icobi_plotter.add_icohps(
            label="NaCl", icohpcollection=self.icobilist_nacl.icohpcollection
        )

        ref_xdata = self.icobilist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icobilist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(
                sum(ydata.values())
            )  # get absolute icohp values as in plots

        fig_xydata = (
            self.icobi_plotter.get_plot()
            .gcf()
            .axes[0]
            .get_children()[0]
            .get_offsets()
            .data
        )

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        self.assertListEqual(ref_xdata, fig_xdata)
        self.assertListEqual(ref_ydata, fig_ydata)

        fig_xydata = self.icobi_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        self.assertListEqual([0, 4], fig_x_lims)
        self.assertListEqual([0, 6], fig_y_lims)

        # icoop

        self.icoop_plotter = IcohpDistancePlotter(are_coops=True)
        self.icoop_plotter.add_icohps(
            label="NaCl", icohpcollection=self.icooplist_nacl.icohpcollection
        )

        ref_xdata = self.icooplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icooplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(
                sum(ydata.values())
            )  # get absolute icohp values as in plots

        fig_xydata = (
            self.icoop_plotter.get_plot()
            .gcf()
            .axes[0]
            .get_children()[0]
            .get_offsets()
            .data
        )

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        self.assertListEqual(ref_xdata, fig_xdata)
        self.assertListEqual(ref_ydata, fig_ydata)

        fig_xydata = self.icoop_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        self.assertListEqual([0, 4], fig_x_lims)
        self.assertListEqual([0, 6], fig_y_lims)


class TestPlotterExceptions(unittest.TestCase):
    def test_plotter_exception(self):
        with self.assertRaises(Exception) as err:
            self.iplotter = InteractiveCohpPlotter()

            data = {"N4: 1 x N-N": []}

            self.iplotter.add_cohps_from_plot_data(plot_data_dict=data, suffix="")

            self.assertEqual(
                err.exception.__str__(),
                "The data provided could not be converted to cohp object.Please recheck the input data",
            )

        with self.assertRaises(Exception) as err:
            self.iplotter = InteractiveCohpPlotter(are_cobis=True, are_coops=True)

            fig = self.iplotter.get_plot()

            self.assertEqual(
                err.exception.__str__(),
                "Plot data should not contain COBI and COOP data at same time",
            )

        with self.assertRaises(Exception) as err:
            self.plotter = PlainCohpPlotter(are_cobis=True, are_coops=True)

            fig = self.plotter.get_plot()

            self.assertEqual(
                err.exception.__str__(),
                "Plot data should not contain COBI and COOP data at same time",
            )
