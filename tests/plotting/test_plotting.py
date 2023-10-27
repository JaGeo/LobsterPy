from __future__ import annotations

import pytest
import gzip
import json
import numpy as np
from pathlib import Path
from plotly.io import read_json
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.io.lobster import Doscar, Icohplist
from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from lobsterpy.plotting import (
    PlainCohpPlotter,
    InteractiveCohpPlotter,
    IcohpDistancePlotter,
    PlainDosPlotter,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestInteractiveCohpPlotter:
    def setup_method(self):
        self.analyse_NaCl = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            summed_spins=False,
        )

        self.analyse_CdF_orb = Analysis(
            path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
            orbital_cutoff=0.10,
            summed_spins=False,
            orbital_resolved=True,
        )

        self.analyse_NaCl_cobi = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="cation-anion",
            cutoff_icohp=0.1,
            summed_spins=False,
            noise_cutoff=0.001,
            are_cobis=True,
        )

        self.analyse_NaCl_cobi_orb = Analysis(
            path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
            summed_spins=False,
            noise_cutoff=0.001,
            orbital_resolved=True,
            are_cobis=True,
        )

        self.analyse_NaSi = Analysis(
            path_to_poscar=TestDir / "test_data/NaSi/POSCAR",
            path_to_cohpcar=TestDir / "test_data/NaSi/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/NaSi/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/NaSi/CHARGE.lobster",
            which_bonds="all",
            cutoff_icohp=0.1,
            summed_spins=True,
        )

        self.analyse_BaTiO3_orb = Analysis(
            path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
            which_bonds="all",
            summed_spins=False,
            orbital_cutoff=0.10,
            orbital_resolved=True,
        )

        self.analyse_K3Sb = Analysis(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            which_bonds="all",
            cutoff_icohp=0.1,
            summed_spins=False,
        )

        plot_data_file_name = TestDir / "test_data/interactive_plotter_ref/mp-8818.json.gz"

        with gzip.open(plot_data_file_name, "rb") as f:
            data = json.loads(f.read().decode("utf-8"))

        lobsterpy_plot_data = {}
        for item in data:
            lobsterpy_plot_data.update(item)

        self.lobsterpy_plot_data = lobsterpy_plot_data["all_bonds"]["lobsterpy_data"]["cohp_plot_data"]

    def test_add_all_relevant_cohps_na_cl(self):
        self.iplotter = InteractiveCohpPlotter(zero_at_efermi=False)

        self.iplotter.add_all_relevant_cohps(analyse=self.analyse_NaCl, label_resolved=False, suffix="")
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 1

        fig = self.iplotter.get_plot(invert_axes=False)
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

    def test_add_all_relevant_cohps_na_cl_cobi_orb(self):
        self.iplotter = InteractiveCohpPlotter(are_cobis=True)

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_NaCl_cobi_orb,
            label_resolved=False,
            suffix="",
            orbital_resolved=True,
        )
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 3

        fig = self.iplotter.get_plot()
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

    def test_add_all_relevant_cohps_ba_ti_o3_orb(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_BaTiO3_orb,
            label_resolved=True,
            orbital_resolved=True,
            suffix="",
        )
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 6

        fig = self.iplotter.get_plot()
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

    def test_add_all_relevant_cohps_cd_f(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(
            analyse=self.analyse_CdF_orb,
            label_resolved=False,
            orbital_resolved=True,
            suffix="",
        )
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 3

        fig = self.iplotter.get_plot()
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

    def test_add_all_relevant_cohps_na_cl_cobi(self):
        self.iplotter = InteractiveCohpPlotter(zero_at_efermi=False, are_cobis=True)

        self.iplotter.add_all_relevant_cohps(analyse=self.analyse_NaCl_cobi, label_resolved=False, suffix="")
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 1

        fig = self.iplotter.get_plot()
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

    def test_add_all_relevant_cohps_k3_sb(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(analyse=self.analyse_K3Sb, label_resolved=True, suffix="")
        assert "All" in self.iplotter._cohps
        assert "K1: 8 x K-K" in self.iplotter._cohps
        assert "K1: 6 x K-Sb" in self.iplotter._cohps
        assert "K2: 4 x K-Sb" in self.iplotter._cohps
        assert "K2: 10 x K-K" in self.iplotter._cohps
        assert "K2: 10 x K-K" in self.iplotter._cohps
        assert "Sb4: 14 x K-Sb" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 6

        fig = self.iplotter.get_plot(sigma=0.3, xlim=[-5, 5], ylim=[-10, 10])
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

    def test_add_cohps_by_lobster_label_na_cl(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_cohps_by_lobster_label(analyse=self.analyse_NaCl, label_list=["5", "10", "15"], suffix="")
        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 1

        fig = self.iplotter.get_plot(integrated=True)
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

    def test_add_cohps_from_plot_data(self):
        self.des = Description(analysis_object=self.analyse_NaSi)

        fig = self.des.plot_interactive_cohps(hide=True)
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

    def test_add_cohps_from_plot_data_json(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_cohps_from_plot_data(plot_data_dict=self.lobsterpy_plot_data, suffix="")

        assert "All" in self.iplotter._cohps
        assert len(self.iplotter._cohps) == 1

        fig = self.iplotter.get_plot()
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

    def test_plot_colors(self):
        self.iplotter = InteractiveCohpPlotter()

        self.iplotter.add_all_relevant_cohps(analyse=self.analyse_K3Sb, label_resolved=True, suffix="")

        Fig_ref = self.iplotter.get_plot()

        color = ["#00FFFF", "#008080", "#00008B", "#808000"]
        # custom color plot
        Fig_cust_col = self.iplotter.get_plot(colors=color)

        for ref, cust in zip(Fig_ref.data, Fig_cust_col.data):
            assert ref.line.color != cust.line.color

    def test_plot_labels(self):
        # plain cohp plotter
        self.plotter = PlainCohpPlotter(are_cobis=True)
        fig = self.plotter.get_plot().gca()

        assert fig.get_xlabel() == "COBI"

        self.plotter = PlainCohpPlotter(are_coops=True)
        fig = self.plotter.get_plot().gca()

        assert fig.get_xlabel() == "COOP"

        self.plotter = PlainCohpPlotter()
        fig = self.plotter.get_plot().gca()

        assert fig.get_xlabel() == "$-$COHP (eV)"

        # interactive plotter
        self.iplotter = InteractiveCohpPlotter(are_cobis=True)
        fig = self.iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "COBI"

        self.iplotter = InteractiveCohpPlotter(are_coops=True)
        fig = self.iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "COOP"

        self.iplotter = InteractiveCohpPlotter()
        fig = self.iplotter.get_plot()

        assert fig.layout.xaxis["title"]["text"] == "−COHP (eV)"
        assert fig.layout.yaxis["title"]["text"] == "$E - E_f \\text{ (eV)}$"

    def test_plaincohp_plotter_options(self):
        self.plotter = PlainCohpPlotter(zero_at_efermi=False)

        for label, cohp in self.lobsterpy_plot_data.items():
            cohp_obj = Cohp.from_dict(cohp)
            self.plotter.add_cohp(label=label, cohp=cohp_obj)

        fig = self.plotter.get_plot(integrated=True, xlim=(-5, 2), ylim=(-4, 4), invert_axes=False).gca()
        assert fig.get_ylabel() == "$-$ICOHP (eV)"
        assert fig.get_xlabel() == "$E$ (eV)"


class TestIcohpDistancePlotter:
    def setup_method(self):
        self.icohplist_nacl = Icohplist(filename=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz")
        self.icobilist_nacl = Icohplist(
            filename=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
            are_cobis=True,
        )
        self.icooplist_nacl = Icohplist(
            filename=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz",
            are_coops=True,
        )

    def test_icohp_plotter_labels(self):
        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(label="NaCl", icohpcollection=self.icohplist_nacl.icohpcollection)
        fig = self.icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "$-$" + "ICOHP (eV)"

        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(label="NaCl_icohp", icohpcollection=self.icohplist_nacl.icohpcollection)
        fig = self.icohp_plotter.get_plot(plot_negative=True).gca()
        assert fig.get_ylabel() == "$-$" + "ICOHP (eV)"

        self.icohp_plotter = IcohpDistancePlotter(are_cobis=True)
        self.icohp_plotter.add_icohps(label="NaCl_icobi", icohpcollection=self.icobilist_nacl.icohpcollection)
        fig = self.icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "ICOBI"

        self.icohp_plotter = IcohpDistancePlotter(are_coops=True)
        self.icohp_plotter.add_icohps(label="NaCl_icoop", icohpcollection=self.icooplist_nacl.icohpcollection)
        fig = self.icohp_plotter.get_plot().gca()
        assert fig.get_ylabel() == "ICOOP"
        assert fig.get_xlabel() == "Bond lengths (Å)"

    def test_plot_data(self):
        self.icohp_plotter = IcohpDistancePlotter()
        self.icohp_plotter.add_icohps(label="NaCl", icohpcollection=self.icohplist_nacl.icohpcollection)

        ref_xdata = self.icohplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icohplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(abs(sum(ydata.values())))  # get absolute icohp values as in plots

        fig_xydata = self.icohp_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = self.icohp_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert [0, 4] == fig_x_lims
        assert [0, 6] == fig_y_lims

        # icobi

        self.icobi_plotter = IcohpDistancePlotter(are_cobis=True)
        self.icobi_plotter.add_icohps(label="NaCl", icohpcollection=self.icobilist_nacl.icohpcollection)

        ref_xdata = self.icobilist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icobilist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(sum(ydata.values()))  # get absolute icohp values as in plots

        fig_xydata = self.icobi_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = self.icobi_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert [0, 4] == fig_x_lims
        assert [0, 6] == fig_y_lims

        # icoop

        self.icoop_plotter = IcohpDistancePlotter(are_coops=True)
        self.icoop_plotter.add_icohps(label="NaCl", icohpcollection=self.icooplist_nacl.icohpcollection)

        ref_xdata = self.icooplist_nacl.icohpcollection._list_length
        ref_ydata = []
        for ydata in self.icooplist_nacl.icohpcollection._list_icohp:
            ref_ydata.append(sum(ydata.values()))  # get absolute icohp values as in plots

        fig_xydata = self.icoop_plotter.get_plot().gcf().axes[0].get_children()[0].get_offsets().data

        fig_xdata = [row[0] for row in fig_xydata]
        fig_ydata = [row[1] for row in fig_xydata]

        assert ref_xdata == fig_xdata
        assert ref_ydata == fig_ydata

        fig_xydata = self.icoop_plotter.get_plot(xlim=(0, 4), ylim=(0, 6)).gcf()

        fig_x_lims = list(fig_xydata.axes[0].get_children()[5].get_view_interval())
        fig_y_lims = list(fig_xydata.axes[0].get_children()[6].get_view_interval())

        assert [0, 4] == fig_x_lims
        assert [0, 6] == fig_y_lims


class TestPlotterExceptions:
    def test_plotter_exception(self):
        with pytest.raises(Exception) as err:
            self.iplotter = InteractiveCohpPlotter()

            data = {"N4: 1 x N-N": []}

            self.iplotter.add_cohps_from_plot_data(plot_data_dict=data, suffix="")

        assert str(err.value) == "The data provided could not be converted to cohp object.Please recheck the input data"

        with pytest.raises(Exception) as err:
            self.iplotter = InteractiveCohpPlotter(are_cobis=True, are_coops=True)

            _ = self.iplotter.get_plot()

        assert str(err.value) == "Plot data should not contain COBI and COOP data at same time"

        with pytest.raises(Exception) as err:
            self.plotter = PlainCohpPlotter(are_cobis=True, are_coops=True)

            _ = self.plotter.get_plot()

        assert str(err.value) == "Plot data should not contain COBI and COOP data at same time"


class TestPlainDosPlotter:
    def setup_method(self) -> None:
        self.NaCl_dos = Doscar(
            doscar=TestDir / "test_data/NaCl_comp_range/DOSCAR.LSO.lobster.gz",
            structure_file=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        )

        self.K3Sb_dos = Doscar(
            doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
            structure_file=TestDir / "test_data/K3Sb/POSCAR.gz",
        )

    def test_na_cl_dos(self):
        complete_dos_obj = self.NaCl_dos.completedos
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

        complete_dos_obj_norm = self.NaCl_dos.completedos.get_normalized()
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

    def test_k3_sb_dos(self):
        complete_dos_obj = self.K3Sb_dos.completedos
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

        # add and test total non normalized smeared dos data and axis labels in the plot
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

    def test_dos_plotter_exceptions(self):
        with pytest.raises(ValueError) as err:
            self.dp = PlainDosPlotter(summed=True, stack=False, sigma=None)

            _ = self.dp.add_site_orbital_dos(site_index=0, orbital="5_s", dos=self.NaCl_dos.completedos)

        assert (
            str(err.value)
            == "Requested orbital is not available for this site, available orbitals are 3s, 2p_y, 2p_z, 2p_x"
        )
