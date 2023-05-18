from __future__ import annotations

import unittest
from pathlib import Path
from plotly.io import read_json, write_json
from lobsterpy.cohp.analyze import Analysis
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

        self.analyse_K3Sb = Analysis(
            path_to_poscar=TestDir / "TestData/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/K3Sb/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/K3Sb/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/K3Sb/CHARGE.lobster.gz",
            whichbonds="all",
            cutoff_icohp=0.1,
            summed_spins=False,
        )

    def test_add_all_relevant_cohps_NaCl(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_NaCl, label_resolved=False, label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertIn("Na1: 6 x Cl-Na", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 3)

        Fig = self.iplotter.get_plot(invert_axes=False)
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_NaCl.json",
            engine="json",
        )

        self.assertEqual(Fig.data, ref_fig.data)

    def test_add_all_relevant_cohps_K3Sb(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_K3Sb, label_resolved=True, label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertIn("K1: 8 x K-K", self.iplotter._cohps)
        self.assertIn("K1: 6 x Sb-K", self.iplotter._cohps)
        self.assertIn("K2: 4 x Sb-K", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("K2: 10 x K-K", self.iplotter._cohps)
        self.assertIn("Sb4: 14 x Sb-K", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 7)

        Fig = self.iplotter.get_plot(sigma=0.3, xlim=[-5, 5], ylim=[-10, 10])
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_K3Sb.json",
            engine="json",
        )

        self.assertEqual(Fig.data, ref_fig.data)

    def test_add_cohps_by_lobster_label_NaCl(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_cohps_by_lobster_label(
            analyse=self.analyse_NaCl, label_list=["5", "10", "15"], label_addition=""
        )
        self.assertIn("Please select COHP label here", self.iplotter._cohps)
        self.assertIn("All", self.iplotter._cohps)
        self.assertIn("Na-Na: 5", self.iplotter._cohps)
        self.assertIn("Na-Na: 10", self.iplotter._cohps)
        self.assertIn("Na-Na: 15", self.iplotter._cohps)
        self.assertEqual(len(self.iplotter._cohps), 5)

        Fig = self.iplotter.get_plot(integrated=True)
        ref_fig = read_json(
            TestDir / "TestData/interactive_plotter_ref/analyse_NaCl_label.json",
            engine="json",
        )

        self.assertEqual(Fig.data, ref_fig.data)
