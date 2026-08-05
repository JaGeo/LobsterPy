import warnings
from pathlib import Path

import pytest
from pymatgen.core import Structure
from pymatgen.io.lobster import Bandoverlaps, Charge, Doscar, Lobsterin, Lobsterout
from pymatgen.io.vasp import Vasprun

from lobsterpy.quality import LobsterCalcQuality

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestLobsterCalcQuality:
    def test_calc_quality_summary_consistency(self):
        charge_obj = Charge(filename=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz")
        bandoverlaps_obj = Bandoverlaps(filename=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz")
        structure_obj = Structure.from_file(filename=TestDir / "test_data" / "K3Sb" / "CONTCAR.gz")
        vasprun_obj = Vasprun(
            filename=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz", parse_eigen=False, parse_potcar_file=False
        )
        doscar = Doscar(
            doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            structure_file=None,
            structure=structure_obj,
        )
        lobsterin_obj = Lobsterin.from_file(TestDir / "test_data" / "K3Sb" / "lobsterin.gz")
        lobsterout_obj = Lobsterout(filename=TestDir / "test_data" / "K3Sb" / "lobsterout.gz")

        calc_quality_with_objs = LobsterCalcQuality(
            structure=structure_obj,
            lobster_dos=doscar.completedos,
            charge=charge_obj,
            potcar_symbols=["K_sv", "Sb"],
            bandoverlaps=bandoverlaps_obj,
            vasp_dos=vasprun_obj.complete_dos,
            lobsterin=lobsterin_obj,
            lobsterout=lobsterout_obj,
        )

        calc_quality_dict_with_objs = calc_quality_with_objs.get_calculation_quality_summary(
            dos_comparison=True, bva_comp=True, n_bins=256
        )

        calc_quality_from_paths = LobsterCalcQuality.from_files(
            poscar=TestDir / "test_data" / "K3Sb" / "CONTCAR.gz",
            potcar_symbols=["K_sv", "Sb"],
            charge=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz",
            doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            vasprun=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz",
            bandoverlaps=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz",
            lobsterout=TestDir / "test_data" / "K3Sb" / "lobsterout.gz",
            lobsterin=TestDir / "test_data" / "K3Sb" / "lobsterin.gz",
        )

        calc_quality_dict_from_paths = calc_quality_from_paths.get_calculation_quality_summary(
            dos_comparison=True, bva_comp=True, n_bins=256
        )

        assert calc_quality_dict_with_objs == calc_quality_dict_from_paths

        # Test new LobsterCalcQuality class method
        calc_quality_dict_from_dir = LobsterCalcQuality.from_directory(
            path_to_lobster_calc=TestDir / "test_data" / "K3Sb"
        ).get_calculation_quality_summary(dos_comparison=True, bva_comp=True, n_bins=256)

        assert calc_quality_dict_with_objs == calc_quality_dict_from_dir

    def test_calc_quality_summary_exceptions(self):
        charge_obj = Charge(filename=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz")
        bandoverlaps_obj = Bandoverlaps(filename=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz")
        structure_obj = Structure.from_file(filename=TestDir / "test_data" / "K3Sb" / "CONTCAR.gz")
        vasprun_obj = Vasprun(
            filename=TestDir / "test_data" / "K3Sb" / "vasprun.xml.gz", parse_eigen=False, parse_potcar_file=False
        )
        doscar = Doscar(
            doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
            structure_file=None,
            structure=structure_obj,
        )
        lobsterin_obj = Lobsterin.from_file(TestDir / "test_data" / "K3Sb" / "lobsterin.gz")
        lobsterout_obj = Lobsterout(filename=TestDir / "test_data" / "K3Sb" / "lobsterout.gz")

        with pytest.raises(ValueError) as err_potcar:  # noqa: PT011
            self.calc_des = LobsterCalcQuality.from_files(
                poscar=TestDir / "test_data" / "K3Sb" / "CONTCAR.gz",
                potcar=None,
                potcar_symbols=None,
                doscar=TestDir / "test_data" / "K3Sb" / "DOSCAR.LSO.lobster.gz",
                charge=TestDir / "test_data" / "K3Sb" / "CHARGE.lobster.gz",
                bandoverlaps=TestDir / "test_data" / "K3Sb" / "bandOverlaps.lobster.gz",
                vasprun=None,
                lobsterin=TestDir / "test_data" / "K3Sb" / "lobsterin.gz",
                lobsterout=TestDir / "test_data" / "K3Sb" / "lobsterout.gz",
            )

        assert (
            str(err_potcar.value) == "Provide potcar, potcar_symbols, or vasprun.xml to "
            "identify basis set used for projections"
        )

        with pytest.raises(Exception) as err_charge:  # noqa: PT011
            self.calc_des = LobsterCalcQuality(
                structure=structure_obj,
                lobster_dos=doscar.completedos,
                charge=None,
                potcar_symbols=["K_sv", "Sb"],
                bandoverlaps=bandoverlaps_obj,
                vasp_dos=vasprun_obj.complete_dos,
                lobsterin=lobsterin_obj,
                lobsterout=lobsterout_obj,
            ).get_calculation_quality_summary(bva_comp=True)

        assert str(err_charge.value) == "BVA comparison requested but CHARGE.lobster not provided"

        with pytest.raises(ValueError) as err_doscar:  # noqa: PT011
            self.calc_des = LobsterCalcQuality(
                structure=structure_obj,
                lobster_dos=None,
                charge=charge_obj,
                potcar_symbols=["K_sv", "Sb"],
                bandoverlaps=bandoverlaps_obj,
                vasp_dos=None,
                lobsterin=lobsterin_obj,
                lobsterout=lobsterout_obj,
            ).get_calculation_quality_summary(dos_comparison=True, n_bins=256)

        assert str(err_doscar.value) == "DOS comparison requested but DOS data missing"


class TestCalcQualityDescribe:
    def test_calc_quality_description_text(self):
        calc_quality_K3Sb = LobsterCalcQuality.from_files(
            poscar=TestDir / "test_data/K3Sb/CONTCAR.gz",
            charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            lobsterout=TestDir / "test_data/K3Sb/lobsterout.gz",
            lobsterin=TestDir / "test_data/K3Sb/lobsterin.gz",
            potcar_symbols=["K_sv", "Sb"],
            bandoverlaps=TestDir / "test_data/K3Sb/bandOverlaps.lobster.gz",
            doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
            vasprun=TestDir / "test_data/K3Sb/vasprun.xml.gz",
        )

        calc_quality_K3Sb_dict = calc_quality_K3Sb.get_calculation_quality_summary(
            dos_comparison=True, bva_comp=True, n_bins=256, e_range=[-20, 0]
        )

        calc_quality_CsH = LobsterCalcQuality.from_files(
            poscar=TestDir / "test_data/CsH/CONTCAR.gz",
            charge=TestDir / "test_data/CsH/CHARGE.lobster.gz",
            lobsterout=TestDir / "test_data/CsH/lobsterout.gz",
            lobsterin=TestDir / "test_data/CsH/lobsterin.gz",
            potcar_symbols=["Cs_sv", "H"],
            bandoverlaps=TestDir / "test_data/CsH/bandOverlaps.lobster.gz",
        )

        calc_quality_CsH_dict = calc_quality_CsH.get_calculation_quality_summary(dos_comparison=False, bva_comp=True)

        calc_quality_k3sb_des_text = LobsterCalcQuality.describe(calc_quality_K3Sb_dict)
        assert calc_quality_k3sb_des_text == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 0.83 and 6.36 %, respectively.",
            (
                "The bandOverlaps.lobster file is generated during the LOBSTER run. This indicates that the "
                "projected wave function is not completely orthonormalized; however, the maximal deviation values "
                "observed compared to the identity matrix is below the threshold of 0.1."
            ),
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis.",
            (
                "The Tanimoto index from DOS comparisons in the energy range between -20, 0 eV for s, p,"
                " summed orbitals are: 0.8532, 0.9481, 0.9275."
            ),
        ]

        calc_quality_csh_des_text = LobsterCalcQuality.describe(calc_quality_CsH_dict)
        assert calc_quality_csh_des_text == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 3.01 and 13.73 %, respectively.",
            (
                "The bandOverlaps.lobster file is generated during the LOBSTER run. This indicates that the "
                "projected wave function is not completely orthonormalized. "
                "The maximal deviation value from the identity matrix is 0.4285, and there are 0.1822 percent "
                "k-points above the deviation threshold of 0.1. Please check the results of other quality checks "
                "like dos comparisons, charges, charge spillings before using the results for further analysis."
            ),
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis.",
        ]


class TestCalcQualityDescribeWarnings:
    def test_warnings(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("once")
            calc_quality_obj = LobsterCalcQuality.from_files(
                poscar=TestDir / "test_data/BaTe_low_quality/POSCAR.lobster.vasp.gz",
                charge=TestDir / "test_data/BaTe_low_quality/CHARGE.lobster.gz",
                lobsterout=TestDir / "test_data/BaTe_low_quality/lobsterout.gz",
                lobsterin=TestDir / "test_data/BaTe_low_quality/lobsterin.gz",
                potcar_symbols=["Ba_sv", "Te"],
                doscar=TestDir / "test_data/BaTe_low_quality/DOSCAR.lobster.gz",
                vasprun=TestDir / "test_data/BaTe_low_quality/vasprun.xml.gz",
                # e_range=[-50, 60],
                # dos_comparison=True,
                # bva_comp=False,
                # n_bins=500,
            )

            calc_quality_warnings = calc_quality_obj.get_calculation_quality_summary(
                dos_comparison=True, bva_comp=False, n_bins=500, e_range=[-50, 60]
            )

        # messages = []
        actual_warnings = [str(warning.message) for warning in w]

        expected_warnings = [
            (
                "This method is being deprecated and will be removed on 30-06-2026. "
                "Please use `lobsterpy.quality.LobsterCalcQuality.from_files()` or "
                "`lobsterpy.quality.LobsterCalcQuality.from_directory()` instead."
            ),
            "Consider using DOSCAR.LSO.lobster, as non LSO DOS from LOBSTER can have negative DOS values",
            (
                "Minimum energy range requested for DOS comparisons is not available in VASP or LOBSTER calculation. "
                "Thus, setting `min_e` to the minimum possible value of -15 eV"
            ),
            (
                "Maximum energy range requested for DOS comparisons is not available in VASP or LOBSTER calculation. "
                "Thus, setting `max_e` to the maximum possible value of 5 eV"
            ),
            (
                "Number of bins requested for DOS comparisons is larger than the number of points in the energy "
                "interval. Thus, setting `n_bins` to 107."
            ),
            (
                "Input DOS files have very few points in the energy interval and thus comparisons will not be "
                "reliable. Please rerun the calculations with higher number of DOS points. "
                "Set NEDOS and COHPSteps tags to >= 2000 in VASP and LOBSTER calculations, respectively."
            ),
        ]

        for actual, expected in zip(actual_warnings, expected_warnings):
            assert actual == expected

        calc_des = LobsterCalcQuality.describe(calc_quality_warnings)

        assert calc_des == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 2.255 and 12.72 %, respectively.",
            (
                "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
                "generated during the LOBSTER run."
            ),
            (
                "The Tanimoto index from DOS comparisons in the energy range between -15, 5 eV for s, p, summed "
                "orbitals are: 0.6712, 0.8113, 0.8064."
            ),
        ]

        with warnings.catch_warnings(record=True) as w2:
            warnings.simplefilter("once")
            calc_quality_warnings2_obj = LobsterCalcQuality.from_files(
                poscar=TestDir / "test_data/C/CONTCAR.gz",
                charge=TestDir / "test_data/C/CHARGE.lobster.gz",
                lobsterout=TestDir / "test_data/C/lobsterout.gz",
                lobsterin=TestDir / "test_data/C/lobsterin.gz",
                potcar_symbols=["C"],
            )
            calc_quality_warnings2 = calc_quality_warnings2_obj.get_calculation_quality_summary(bva_comp=True)

        assert "Oxidation states from BVA analyzer cannot" in str(w2[-1].message)

        calc_des2 = LobsterCalcQuality.describe(calc_quality_warnings2)

        assert calc_des2 == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 0.98 and 8.93 %, respectively.",
            (
                "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
                "generated during the LOBSTER run."
            ),
            "Oxidation states from BVA analyzer cannot be determined. Thus BVA charge comparison is not conducted.",
        ]

        with warnings.catch_warnings(record=True) as w3:
            warnings.simplefilter("once")
            warnings.filterwarnings("ignore", module="pymatgen")
            warnings.filterwarnings("ignore", module="spglib")
            calc_quality_warnings3_obj = LobsterCalcQuality.from_files(
                poscar=TestDir / "test_data/BeTe/CONTCAR.gz",
                charge=TestDir / "test_data/BeTe/CHARGE.lobster.gz",
                lobsterout=TestDir / "test_data/BeTe/lobsterout.gz",
                lobsterin=TestDir / "test_data/BeTe/lobsterin.gz",
                potcar_symbols=["Be_sv", "Te"],
            )

            calc_quality_warnings3 = calc_quality_warnings3_obj.get_calculation_quality_summary(bva_comp=True)
        assert (
            str(w3[0].message) == "Consider rerunning the calc with the minimum basis as well. "
            "Choosing a larger basis set is recommended if charge "
            "spilling significantly improves."
        )
        calc_des3 = LobsterCalcQuality.describe(calc_quality_warnings3)

        assert calc_des3 == [
            (
                "Consider rerunning the calculation with the minimum basis as well. "
                "Choosing a larger basis set is only recommended if you see a significant improvement of "
                "the charge spilling."
            ),
            "The absolute and total charge spilling for the calculation is 1.48 and 13.99 %, respectively.",
            (
                "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
                "generated during the LOBSTER run."
            ),
            "The atomic charge signs from Mulliken population analysis do not agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis do not agree with the bond valence analysis.",
        ]

        with warnings.catch_warnings(record=True) as w5:
            warnings.simplefilter("once")
            calc_quality_warnings5 = LobsterCalcQuality.from_files(
                poscar=TestDir / "test_data/C/CONTCAR.gz",
                charge=TestDir / "test_data/C/CHARGE.lobster.gz",
                lobsterout=TestDir / "test_data/C/lobsterout.gz",
                lobsterin=TestDir / "test_data/C/lobsterin.gz",
                potcar_symbols=["C"],
            )
            _summary_dict = calc_quality_warnings5.get_calculation_quality_summary(bva_comp=True)
        assert "Oxidation states from BVA analyzer cannot" in str(w5[-1].message)
