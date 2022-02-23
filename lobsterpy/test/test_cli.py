import unittest
import os
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis
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

    def tearDown(self) -> None:
        os.chdir(CurrentDir)
