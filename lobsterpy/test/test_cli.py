from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
from typing import Union

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.style
import pytest

from lobsterpy.cli import get_parser, run

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"
ref_data_file = TestDir / "TestData/cli-reference.json"
test_cases = [
    ["automatic-plot"],
    ["automaticplot", "--allbonds"],
    ["automatic-plot", "--style", "ggplot"],
    ["description"],
    ["description", "--all-bonds"],
    ["plot", "1", "2"],
    ["plot", "1", "--cobis"],
    ["plot", "1", "--coops"],
    ["plot", "--summed", "1", "2"],
    ["plot", "1", "2", "--orbitalwise", "3s-3s", "3s-3s"],
    ["plot", "1", "--orbitalwise", "all"],
    ["plot", "1", "--fontsize", "20"],
    ["plot", "1", "--width", "20", "--height", "20"],
    ["plot", "1", "--width", "20"],
    ["plot", "1", "--height", "20"],
    ["plot", "1", "--style", "dark_background"],
]

error_test_cases = [
    (["plot", "1", "2", "--orbitalwise", "3s-3s"], IndexError),
    (["plot", "400000"], IndexError),
    (["plot", "1", "--orbitalwise", "1s-1s"], IndexError),
]


class TestCLI:
    @pytest.fixture
    def inject_mocks(self, mocker):
        # Disable calls to show() so we can get the current figure using gcf()
        mocker.patch("matplotlib.pyplot.show")
        mocker.resetall()

    @pytest.fixture
    def clean_plot(self):
        yield
        plt.close("all")
        matplotlib.style.use("default")

    with open(ref_data_file, "r") as fd:
        ref_results = json.load(fd)

    @classmethod
    def setup_class(cls):
        os.chdir(TestDir / "TestData/NaCl")

    @pytest.mark.parametrize("args", test_cases)
    def test_cli_results(self, args, capsys, inject_mocks, clean_plot):
        test = get_parser().parse_args(args)
        run(test)

        captured = capsys.readouterr()
        assert captured.out == self.ref_results[" ".join(args)]["stdout"]

        plot_attributes = self.get_plot_attributes(plt.gcf())
        ref_plot_attributes = self.ref_results[" ".join(args)]["plot"]

        if ref_plot_attributes is None:
            assert plot_attributes is None

        else:
            for key, ref_value in ref_plot_attributes.items():
                if key == "xydata":
                    for line, ref_line in zip(plot_attributes[key], ref_value):
                        assert np.array(np.array(line)) == pytest.approx(
                            np.array(ref_line)
                        )
                else:
                    assert plot_attributes[key] == pytest.approx(ref_value)

    @pytest.mark.parametrize("args, error", error_test_cases)
    def test_cli_errors(self, args, error, inject_mocks):
        with pytest.raises(error):
            test = get_parser().parse_args(args)
            run(test)

    def test_plot_saved(self, tmp_path, inject_mocks, clean_plot):
        plot_path = tmp_path / "plot.png"
        args = ["plot", "1", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

    def test_json_saved(self, tmp_path, inject_mocks, clean_plot):
        json_path = tmp_path / "data.json"
        args = ["automatic-plot", "--json", str(json_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(json_path)

    @staticmethod
    def assert_is_finite_file(path: Path) -> None:
        assert path.is_file()
        assert path.stat().st_size > 0

    @pytest.mark.skip(reason="Only enable this test to regenerate test data")
    def test_generate_ref_data(self, inject_mocks):
        json_data = {}

        for args in test_cases:
            plt.close("all")
            matplotlib.style.use("default")

            with redirect_stdout(io.StringIO()) as stdout:
                test = get_parser().parse_args(args)
                run(test)
                json_data[" ".join(args)] = {"stdout": stdout.getvalue()}

                fig = plt.gcf()
                json_data[" ".join(args)].update(
                    {"plot": self.get_plot_attributes(fig)}
                )

        with open(ref_data_file, "w") as fd:
            json.dump(json_data, fd, indent=4, sort_keys=True)

    @staticmethod
    def get_plot_attributes(fig: Figure) -> Union[dict, None]:
        if fig.axes:
            ax = fig.gca()

            return {
                "xydata": [line.get_xydata().tolist() for line in ax.lines],
                "facecolor": ax.get_facecolor(),
                "size": fig.get_size_inches().tolist(),
            }
        else:
            return None

    def teardown_class(self) -> None:
        os.chdir(CurrentDir)
