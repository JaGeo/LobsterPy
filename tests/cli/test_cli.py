from __future__ import annotations

import gzip
import io
import json
import os
import shutil
import sys
from contextlib import redirect_stdout
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.style
import numpy as np
import pytest
from matplotlib.figure import Figure

from lobsterpy.cli import get_parser, run

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"
ref_data_file = TestDir / "test_data/cli-reference.json"
test_cases = [
    ["automatic-plot"],
    ["automaticplot", "--allbonds"],
    ["automatic-plot", "--style", "ggplot"],
    ["automatic-plot", "--sigma", "1.2"],
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
    ["plot", "1", "--sigma", "1.2"],
    ["plot", "1", "--fwhm", "1"],
]

error_test_cases = [
    (["plot", "1", "2", "--orbitalwise", "3s-3s"], IndexError),
    (["plot", "400000"], IndexError),
    (["plot", "1", "--orbitalwise", "1s-1s"], IndexError),
]


class TestCLI:
    @pytest.fixture()
    def inject_mocks(self, mocker):  # noqa: PT004
        # Disable calls to show() so we can get the current figure using gcf()
        mocker.patch("matplotlib.pyplot.show")
        mocker.resetall()

    @pytest.fixture()
    def clean_plot(self):  # noqa: PT004
        yield
        plt.close("all")
        matplotlib.style.use("default")

    with open(ref_data_file) as fd:
        ref_results = json.load(fd)

    @classmethod
    def setup_class(cls):
        os.chdir(TestDir / "test_data/NaCl")

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

    @pytest.mark.parametrize("args, error", error_test_cases)  # noqa: PT006
    def test_cli_errors(self, args, error, inject_mocks):
        with pytest.raises(error):  # noqa: PT012
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
        args = ["automatic-plot", "--file-json", str(json_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(json_path)

    def test_hideplot_cli(self, tmp_path, inject_mocks, clean_plot):
        os.chdir(TestDir / "test_data/NaCl")
        # tests skip showing plots generated using automaticplot
        args = [
            "automaticplot",
            "--hideplot",
        ]
        test = get_parser().parse_args(args)
        run(test)

        # tests tests skip showing plot and directly saves the generated figures
        plot_path = tmp_path / "autoplot.png"
        args = ["automaticplot", "--hideplot", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        # tests plot hide and save case based on label
        plot_path = tmp_path / "plot2.png"
        args = ["plot", "2", "--hideplot", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        # tests plot show and save case based on label
        plot_path = tmp_path / "plot34.png"
        args = ["plot", "3", "4", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

    def test_cli_interactive_plotter(self):
        os.chdir(TestDir / "test_data/NaCl")
        # tests skip showing plots generated using automatic interactive plotter
        args = [
            "automaticplotia",
            "--hideplot",
        ]
        test = get_parser().parse_args(args)
        run(test)

    def test_iaplot_saved(self, tmp_path, inject_mocks, clean_plot):
        plot_path = tmp_path / "plot.html"
        args = ["automaticplotia", "--hideplot", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

    def test_cli_interactive_plotter_cobi(self):
        os.chdir(TestDir / "test_data/NaCl_comp_range")
        # tests skip showing plots generated using automatic interactive plotter
        args = ["automatic-plot-ia", "--hideplot", "--cobis"]
        test = get_parser().parse_args(args)
        run(test)

    def test_cli_interactive_plotter_coops(self):
        os.chdir(TestDir / "test_data/CdF_comp_range")
        # tests skip showing plots generated using automatic interactive plotter
        args = [
            "auto-plot-ia",
            "--orbitalresolved",
            "--hideplot",
            "--coops",
            "--allbonds",
            "--orbitalplot",
        ]
        test = get_parser().parse_args(args)
        run(test)

    def test_icoxxlist_plots(self, tmp_path, inject_mocks, clean_plot):
        os.chdir(TestDir / "test_data/NaCl")
        plot_path = tmp_path / "plot.png"
        args = ["ploticohpdistance", "--hideplot", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        os.chdir(TestDir / "test_data/CdF_comp_range")
        plot_path = tmp_path / "plot.png"
        args = ["ploticohpdistance", "--coops", "--saveplot", str(plot_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        os.chdir(TestDir / "test_data/CdF_comp_range")
        args = ["ploticohpdistance", "--cobis"]
        test = get_parser().parse_args(args)
        run(test)

    def test_lobsterin_generation(self, tmp_path):
        os.chdir(TestDir / "test_data/Test_Input_Generation_Empty")
        lobsterinpath = tmp_path / "lobsterin.lobsterpy"
        INCARpath = tmp_path / "INCAR.lobsterpy"
        args = [
            "create-inputs",
            "-flobsterin",
            str(lobsterinpath),
            "-fincarout",
            str(INCARpath),
        ]
        test = get_parser().parse_args(args)
        run(test)
        for filepath in [
            tmp_path / "lobsterin.lobsterpy-0",
            tmp_path / "lobsterin.lobsterpy-1",
            tmp_path / "INCAR.lobsterpy-0",
            tmp_path / "INCAR.lobsterpy-1",
        ]:
            self.assert_is_finite_file(filepath)

        # test create-inputs alias and overwrite

        os.chdir(TestDir / "test_data/Test_Input_Generation_Empty")
        lobsterinpath = tmp_path / "lobsterin.lobsterpy"
        INCARpath = tmp_path / "INCAR.lobsterpy"
        args = [
            "createinputs",
            "--file-lobsterin",
            str(lobsterinpath),
            "--file-incar-out",
            str(INCARpath),
            "--overwrite",
        ]
        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl")

    def test_lobsterin_generation_error(self, tmp_path):
        os.chdir(TestDir / "test_data/Test_Input_Generation_Empty")
        lobsterinpath = tmp_path / "lobsterin.lobsterpy"
        INCARpath = tmp_path / "INCAR.lobsterpy"
        args = [
            "create-inputs",
            "--file-lobsterin",
            str(lobsterinpath),
            "--file-incar-out",
            str(INCARpath),
        ]
        test = get_parser().parse_args(args)
        run(test)
        for filepath in [
            tmp_path / "lobsterin.lobsterpy-0",
            tmp_path / "lobsterin.lobsterpy-1",
            tmp_path / "INCAR.lobsterpy-0",
            tmp_path / "INCAR.lobsterpy-1",
        ]:
            self.assert_is_finite_file(filepath)

        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            args = [
                "create-inputs",
                "--file-lobsterin",
                str(lobsterinpath),
                "--file-incar-out",
                str(INCARpath),
            ]
            test = get_parser().parse_args(args)
            run(test)
            assert (
                str(err.value)
                == 'please use "--overwrite" if you would like to overwrite existing lobster inputs'
            )

        with pytest.raises(ValueError) as err2:  # noqa: PT012, PT011
            args = [
                "create-inputs",
                "--file-lobsterin",
                str(lobsterinpath),
                "--file-incar-out",
                str(INCARpath),
                "--userbasis",
                "Cd.4d.5s F.2p.2s",
            ]
            test = get_parser().parse_args(args)
            run(test)
            assert (
                str(err2.value)
                == 'please use "--overwrite" if you would like to overwrite existing lobster inputs'
            )

    def test_cli_automatic_analysis_error(self):
        with pytest.raises(ValueError) as err1:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/NaCl")
            args1 = [
                "description",
                "--cobis",
            ]
            test1 = get_parser().parse_args(args1)

            run(test1)
            assert (
                str(err1.value)
                == "Files required for automatic analysis of COBIs (ICOBILIST.lobster and "
                "COBICAR.lobster) not found in the directory"
            )

        with pytest.raises(ValueError) as err2:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/NaCl")
            args2 = [
                "description",
                "--coops",
            ]
            test2 = get_parser().parse_args(args2)
            run(test2)
            assert (
                str(err2.value)
                == "Files required for automatic analysis of COOPs (ICOOPLIST.lobster and "
                "COOPCAR.lobster) not found in the directory"
            )

    def test_lobsterin_generation_error_userbasis(self, tmp_path):
        # This is a test for the user-defined basis set.
        os.chdir(TestDir / "test_data/Test_Input_Generation_Empty")
        lobsterinpath = tmp_path / "lobsterin.lobsterpy"
        INCARpath = tmp_path / "INCAR.lobsterpy"
        args = [
            "create-inputs",
            "--file-lobsterin",
            str(lobsterinpath),
            "--file-incar-out",
            str(INCARpath),
            "--userbasis",
            "Na.3s.3p Cl.3s.3p",
        ]
        test = get_parser().parse_args(args)
        run(test)

        for basis in ["Na 3s 3p", "Cl 3s 3p"]:
            with open(tmp_path / "lobsterin.lobsterpy-0") as f:
                file_data = f.read()
            assert basis in file_data

        for filepath in [
            tmp_path / "lobsterin.lobsterpy-0",
            tmp_path / "INCAR.lobsterpy-0",
        ]:
            self.assert_is_finite_file(filepath)

        os.chdir(TestDir / "test_data/NaCl")

    def test_calc_quality_summary_nacl(self, tmp_path):
        os.chdir(TestDir / "test_data/NaCl_comp_range")
        calc_quality_json_path = tmp_path / "calc_quality_json.json"
        args = [
            "description-quality",
            "--potcar-symbols",
            "Na_pv Cl",
            "--bvacomp",
            "--doscomp",
            "--erange",
            "-20",
            "0",
            "--file-calc-quality-json",
            str(calc_quality_json_path),
        ]
        captured_output = io.StringIO()
        sys.stdout = captured_output

        test = get_parser().parse_args(args)
        run(test)

        calc_quality_text = captured_output.getvalue().strip()

        sys.stdout = sys.__stdout__

        ref_text = (
            "The LOBSTER calculation used minimal basis. "
            "The absolute and total charge spilling for the calculation is 0.3 and 5.58 %, respectively. "
            "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
            "generated during the LOBSTER run. "
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis. "
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis. "
            "The Tanimoto index from DOS comparisons in the energy range between -20, 0 eV for s, p, summed orbitals "
            "are: 0.9966, 0.9977, 0.9822."
        )

        assert calc_quality_text == ref_text
        self.assert_is_finite_file(calc_quality_json_path)

    def test_calc_quality_summary_k3sb(self, tmp_path):
        os.chdir(TestDir / "test_data/K3Sb")
        calc_quality_json_path = tmp_path / "calc_quality_json.json"
        args = [
            "description-quality",
            "--bvacomp",
            "--potcar-symbols",
            "K_sv Sb",
            "--doscomp",
            "--file-doscar",
            "DOSCAR.LSO.lobster",
            "--erange",
            "-20",
            "0",
            "--file-calc-quality-json",
            str(calc_quality_json_path),
        ]
        captured_output = io.StringIO()
        sys.stdout = captured_output

        test = get_parser().parse_args(args)
        run(test)

        calc_quality_text = captured_output.getvalue().strip()

        sys.stdout = sys.__stdout__

        ref_text = (
            "The LOBSTER calculation used minimal basis. "
            "The absolute and total charge spilling for the calculation is 0.83 and 6.36 %, respectively. "
            "The bandOverlaps.lobster file is generated during the LOBSTER run. This indicates that "
            "the projected wave function is not completely orthonormalized; however, the "
            "maximal deviation values observed compared to the identity matrix is below the threshold of 0.1. "
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis. "
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis. "
            "The Tanimoto index from DOS comparisons in the energy range between -20, 0 eV for s, p, summed orbitals "
            "are: 0.8367, 0.9565, 0.9357."
        )

        assert calc_quality_text == ref_text
        self.assert_is_finite_file(calc_quality_json_path)

    def test_dos_plot(self, tmp_path, inject_mocks, clean_plot):
        os.chdir(TestDir / "test_data/K3Sb")
        plot_path = tmp_path / "autoplot.png"

        args = [
            "plot-dos",
            "-addtdos",
            "--spddos",
            "-sspin",
            "--file-doscar",
            "DOSCAR.LSO.lobster",
            "--sigma",
            "0.2",
            "--xlim",
            "-5",
            "0.5",
            "--saveplot",
            str(plot_path),
        ]

        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        plot_path = tmp_path / "autoplot.png"
        args = [
            "plot-dos",
            "--elementdos",
            "--file-doscar",
            "DOSCAR.LSO.lobster",
            "--ylim",
            "-5",
            "0.5",
            "--hideplot",
            "--saveplot",
            str(plot_path),
        ]

        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(plot_path)

        os.chdir(TestDir / "test_data/K3Sb")
        args = [
            "plot-dos",
            "--file-doscar",
            "DOSCAR.LSO.lobster",
            "--element",
            "K",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
            "--site",
            "1",
            "--orbital",
            "3s",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
            "--site",
            "1",
            "--orbital",
            "all",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
            "--site",
            "0",
            "1",
            "--orbital",
            "all",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
            "--site",
            "0",
            "1",
            "--orbital",
            "all",
            "3s",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
            "--site",
            "0",
            "--orbital",
            "all",
            "3s",
            "--invertaxis",
        ]

        test = get_parser().parse_args(args)
        run(test)

        os.chdir(TestDir / "test_data/NaCl_comp_range")
        args = [
            "plot-dos",
        ]

        test = get_parser().parse_args(args)
        run(test)

    def test_cli_exceptions(self):
        # Calc files missing exception test
        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir)
            args = [
                "description",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "Files necessary for automatic analysis of LOBSTER outputs not found in the current directory"
            )

        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir)
            args = [
                "description-quality",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "Mandatory files necessary for LOBSTER calc quality not found in the current directory."
            )

        # doscar comparison exceptions test
        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/NaCl")
            args = [
                "description-quality",
                "--doscomp",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "DOS comparisons requested but DOSCAR.lobster, vasprun.xml file not found."
            )

        # BVA comparison exceptions test
        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/NaCl")
            args = [
                "description-quality",
                "--bvacomp",
                "--file-charge",
                "../CHARGE.lobster",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "BVA charge requested but CHARGE.lobster file not found."
            )

        # Create-inputs exceptions test
        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/CsH")
            args = [
                "create-inputs",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "Files necessary for creating puts for LOBSTER calcs not found in the current directory."
            )

        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/CsH")
            args = [
                "plot-dos",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "DOSCAR.lobster necessary for plotting DOS not found in the current directory."
            )

        with pytest.raises(ValueError) as err:  # noqa: PT012, PT011
            os.chdir(TestDir / "test_data/K3Sb")
            args = [
                "plot-dos",
                "--file-doscar",
                "DOSCAR.LSO.lobster",
                "--site",
                "1",
            ]

            test = get_parser().parse_args(args)
            run(test)

            assert (
                str(err.value)
                == "Please set both args i.e site and orbital to generate the plot"
            )

    def test_nongz_file_cli(self, tmp_path, inject_mocks, clean_plot):
        # test description from gz input files
        os.chdir(TestDir / "test_data/CsH")
        for file in os.listdir():
            shutil.copy(file, tmp_path)
        os.chdir(tmp_path)
        for file in os.listdir():
            if file.endswith(".gz"):
                uncompressed_file_path = file.split(".gz")[0]  # Remove '.gz' extension
                with gzip.open(file, "rb") as f_in, open(
                    uncompressed_file_path, "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    # Delete the source gzipped file
                    os.remove(file)
        args = ["description", "--allbonds"]

        test = get_parser().parse_args(args)
        run(test)

        # test autoplot and json generation from gz input files
        json_path = tmp_path / "data.json"
        args = ["automatic-plot", "--all-bonds", "--file-json", str(json_path)]
        test = get_parser().parse_args(args)
        run(test)
        self.assert_is_finite_file(json_path)

    def test_gz_file_cli_lobsterinput_generation(self, tmp_path):
        os.chdir(TestDir / "test_data/Test_Input_Generation_Empty/gz/")
        lobsterinpath = tmp_path / "lobsterin.lobsterpy"
        INCARpath = tmp_path / "INCAR.lobsterpy"
        args = [
            "create-inputs",
            "--file-lobsterin",
            str(lobsterinpath),
            "--file-incar-out",
            str(INCARpath),
            "--userbasis",
            "Na.3s.3p Cl.3s.3p",
        ]

        test = get_parser().parse_args(args)
        run(test)

        for basis in ["Na 3s 3p", "Cl 3s 3p"]:
            with open(tmp_path / "lobsterin.lobsterpy-0") as f:
                file_data = f.read()
            assert basis in file_data

        for filepath in [
            tmp_path / "lobsterin.lobsterpy-0",
            tmp_path / "INCAR.lobsterpy-0",
        ]:
            self.assert_is_finite_file(filepath)

        os.chdir(TestDir / "test_data/NaCl")

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
    def get_plot_attributes(fig: Figure) -> dict | None:
        if fig.axes:
            ax = fig.gca()

            return {
                "xydata": [line.get_xydata().tolist() for line in ax.lines],  # type: ignore
                "facecolor": ax.get_facecolor(),
                "size": fig.get_size_inches().tolist(),
            }
        else:  # noqa: RET505
            return None

    def teardown_class(self) -> None:
        os.chdir(CurrentDir)
