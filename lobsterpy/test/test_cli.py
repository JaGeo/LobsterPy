import os
import unittest
from pathlib import Path

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
        plot = run(test)

    def test_plot(self):
        test = get_parser().parse_args(["plot", "1", "2"])
        plot = run(test)

    def test_plot_summed(self):
        test = get_parser().parse_args(["plot", "--summed", "1", "2", ])
        plot = run(test)

    def test_plot_orbitalwise(self):
        test = get_parser().parse_args(["plot", "1", "2", "--orbitalwise", "3s-3s", "3s-3s"])
        plot = run(test)

    def test_plot_errors(self):
        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "1", "2", "--orbitalwise", "3s-3s"])
            plot = run(test)

        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "400000"])
            plot = run(test)

        with self.assertRaises(IndexError):
            test = get_parser().parse_args(["plot", "1", "--orbitalwise", "1s-1s"])
            plot = run(test)



    def tearDown(self) -> None:
        os.chdir(CurrentDir)
