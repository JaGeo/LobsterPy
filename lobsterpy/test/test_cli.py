from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path

import pytest

from lobsterpy.cli import get_parser, run

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"
ref_data_file = TestDir / "TestData/cli-reference.json"
test_args = [
    ["automatic-plot"],
    ["automaticplot", "--allbonds"],
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
]

error_test_cases = [
    (["plot", "1", "2", "--orbitalwise", "3s-3s"], IndexError),
    (["plot", "400000"], IndexError),
    (["plot", "1", "--orbitalwise", "1s-1s"], IndexError),
]


# TODO: check if plots are really generated
class TestCLI:
    @pytest.fixture
    def inject_mocks(self, mocker):
        # Prevent calls to show so we can get the current figure using gcf()
        mocker.patch("matplotlib.pyplot.show")
        mocker.resetall()

    with open(ref_data_file, "r") as fd:
        ref_results = json.load(fd)

    @classmethod
    def setup_class(cls):
        os.chdir(TestDir / "TestData/NaCl")

    @pytest.mark.parametrize("args", test_args)
    def test_cli_results(self, args, capsys, inject_mocks):
        test = get_parser().parse_args(args)
        run(test)

        captured = capsys.readouterr()
        assert captured.out == self.ref_results[" ".join(args)]["stdout"]

    @pytest.mark.parametrize("args, error", error_test_cases)
    def test_cli_errors(self, args, error, inject_mocks):
        with pytest.raises(error):
            test = get_parser().parse_args(args)
            run(test)

    @pytest.mark.skip(reason="Only enable this test to regenerate test data")
    def test_generate_ref_data(self, inject_mocks):
        json_data = {}

        for args in test_args:
            with redirect_stdout(io.StringIO()) as stdout:
                test = get_parser().parse_args(args)
                run(test)
                json_data[" ".join(args)] = {"stdout": stdout.getvalue()}

        with open(ref_data_file, "w") as fd:
            json.dump(json_data, fd, indent=4, sort_keys=True)

    # def test_plot_save(self):
    #     with TemporaryDirectory() as d:
    #         temp_file_name = os.path.join(d, "test.pdf")
    #         test = get_parser().parse_args(["plot", "1", "--saveplot", temp_file_name])
    #         run(test)

    # def test_json(self):
    #     with TemporaryDirectory():
    #         test = get_parser().parse_args(["automatic-plot", "--json"])
    #         run(test)

    def teardown_class(self) -> None:
        os.chdir(CurrentDir)
