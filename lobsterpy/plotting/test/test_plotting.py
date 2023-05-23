from __future__ import annotations

import unittest
import gzip
import json
from pathlib import Path
from plotly.io import read_json, write_json
from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from lobsterpy.plotting import PlainCohpPlotter, InteractiveCohpPlotter

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
            analyse=self.analyse_NaCl, label_resolved=False, label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 2)

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
            analyse=self.analyse_K3Sb, label_resolved=True, label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertIn("K1: 8 x K-K", self.iplotter._cohps)
        self.assertIn("K1: 6 x K-Sb", self.iplotter._cohps)
        self.assertIn("K2: 4 x K-Sb", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("Sb4: 14 x K-Sb", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 7)

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
            analyse=self.analyse_NaCl, label_list=["5", "10", "15"], label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 2)

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

        fig = self.des.plot_interactive_cohps(skip_show=True)
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
            plot_data_dict=self.lobsterpy_plot_data, label_addition=""
        )

        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 2)

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


class TestPlotterExceptions(unittest.TestCase):
    def test_plotter_exception(self):
        with self.assertRaises(Exception) as err:
            self.iplotter = InteractiveCohpPlotter()

            data = {"N4: 1 x N-N": []}

            self.iplotter.add_cohps_from_plot_data(
                plot_data_dict=data, label_addition=""
            )

            self.assertEqual(
                err.exception.__str__(),
                "The data provided could not be converted to cohp object.Please recheck the input data",
            )
