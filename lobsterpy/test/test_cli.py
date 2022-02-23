import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from lobsterpy.cli import get_parser, run

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


# TODO: check if text outputs and plots are really generated


class TestDescribe(unittest.TestCase):
    def setUp(self):
        os.chdir(TestDir / "TestData/NaCl")

    def test_automaticplot(self):
        test = get_parser().parse_args(["automatic-plot"])
        run(test)

    def test_automaticplot_allbonds(self):
        test = get_parser().parse_args(["automaticplot", "--allbonds"])
        run(test)

    def test_description(self):
        test = get_parser().parse_args(["description"])
        run(test)

    def test_description_allbonds(self):
        test = get_parser().parse_args(["description", "--allbonds"])
        run(test)

    def test_plot(self):
        test = get_parser().parse_args(["plot", "1", "2"])
        run(test)

    def test_plot_cobis(self):
        test = get_parser().parse_args(["plot", "1", "--cobis"])
        run(test)

    def test_plot_coops(self):
        test = get_parser().parse_args(["plot", "1", "--coops"])
        run(test)

    def test_plot_summed(self):
        test = get_parser().parse_args(
            [
                "plot",
                "--summed",
                "1",
                "2",
            ]
        )
        run(test)

    def test_plot_orbitalwise(self):
        test = get_parser().parse_args(
            ["plot", "1", "2", "--orbitalwise", "3s-3s", "3s-3s"]
        )
        run(test)

    def test_plot_orbitalwise_all(self):
        test = get_parser().parse_args(["plot", "1", "--orbitalwise", "all"])
        run(test)

    def test_plot_errors(self):
        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "1", "2", "--orbitalwise", "3s-3s"])
            run(test)

        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "400000"])
            run(test)

        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "1", "--orbitalwise", "1s-1s"])
            run(test)

    def test_plot_fontsize(self):
        test = get_parser().parse_args(["plot", "1", "--fontsize", "20"])
        run(test)

    def test_plot_width_height(self):
        test = get_parser().parse_args(["plot", "1", "--width", "20", "--height", "20"])
        run(test)

    def test_plot_width(self):
        test = get_parser().parse_args(["plot", "1", "--width", "20"])
        run(test)

    def test_plot_height(self):
        test = get_parser().parse_args(["plot", "1", "--height", "20"])
        run(test)

    def test_plot_save(self):
        with TemporaryDirectory() as d:
            temp_file_name = os.path.join(d, "test.pdf")
            test = get_parser().parse_args(["plot", "1", "--saveplot", temp_file_name])
            run(test)

    def test_json(self):
        with TemporaryDirectory() as d:
            test = get_parser().parse_args(["automatic-plot", "--json"])
            run(test)

    def tearDown(self) -> None:
        os.chdir(CurrentDir)
