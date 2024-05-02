from __future__ import annotations

import warnings
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


class TestDescribe:
    def test_coordination_environment_to_text(self):
        results_dict = {
            "S:1": "single (CN=1)",
            "L:2": "linear (CN=2)",
            "A:2": "angular (CN=2)",
            "TL:3": "trigonal planar (CN=3)",
            "TY:3": "triangular non-coplanar (CN=3)",
            "TS:3": "t-shaped (CN=3)",
            "T:4": "tetrahedral (CN=4)",
            "S:4": "square planar (CN=4)",
            "SY:4": "square non-coplanar (CN=4)",
            "SS:4": "see-saw like (CN=4)",
            "PP:5": "pentagonal (CN=5)",
            "S:5": "square pyramidal (CN=5)",
            "T:5": "trigonal bipyramidal (CN=5)",
            "O:6": "octahedral (CN=6)",
            "T:6": "trigonal prismatic (CN=6)",
            "PP:6": "pentagonal pyramidal (CN=6)",
            "PB:7": "pentagonal bipyramidal (CN=7)",
            "ST:7": "square-face capped trigonal prismatic (CN=7)",
            "ET:7": "end-trigonal-face capped trigonal prismatic (CN=7)",
            "FO:7": "face-capped octahedron (CN=7)",
            "C:8": "cubic (CN=8)",
            "SA:8": "square antiprismatic (CN=8)",
            "SBT:8": "square-face bicapped trigonal prismatic (CN=8)",
            "TBT:8": "triangular-face bicapped trigonal prismatic (CN=8)",
            "DD:8": "dodecahedronal (with triangular faces) (CN=8)",
            "DDPN:8": "dodecahedronal (with triangular faces - p2345 plane normalized) (CN=8)",
            "HB:8": "hexagonal bipyramidal (CN=8)",
            "BO_1:8": "bicapped octahedral (opposed cap faces) (CN=8)",
            "BO_2:8": "bicapped octahedral (cap faces with one atom in common) (CN=8)",
            "BO_3:8": "bicapped octahedral (cap faces with one edge in common) (CN=8)",
            "TC:9": "triangular cupola (CN=9)",
            "TT_1:9": "Tricapped triangular prismatic (three square - face caps) (CN=9)",
            "TT_2:9": "Tricapped triangular prismatic (two square - face caps and one triangular - face cap) (CN=9)",
            "TT_3:9": "Tricapped triangular prism (one square - face cap and two triangular - face caps) (CN=9)",
            "HD:9": "Heptagonal dipyramidal (CN=9)",
            "TI:9": "tridiminished icosohedral (CN=9)",
            "SMA:9": "Square-face monocapped antiprism (CN=9)",
            "SS:9": "Square-face capped square prismatic (CN=9)",
            "TO_1:9": "Tricapped octahedral (all 3 cap faces share one atom) (CN=9)",
            "TO_2:9": "Tricapped octahedral (cap faces are aligned) (CN=9)",
            "TO_3:9": "Tricapped octahedron (all 3 cap faces are sharing one edge of a face) (CN=9)",
            "PP:10": "Pentagonal prismatic (CN=10)",
            "PA:10": "Pentagonal antiprismatic (CN=10)",
            "SBSA:10": "Square-face bicapped square antiprismatic (CN=10)",
            "MI:10": "Metabidiminished icosahedral (CN=10)",
            "S:10": "sphenocoronal (CN=10)",
            "H:10": "Hexadecahedral (CN=10)",
            "BS_1:10": "Bicapped square prismatic (opposite faces) (CN=10)",
            "BS_2:10": "Bicapped square prism(adjacent faces) (CN=10)",
            "TBSA:10": "Trigonal-face bicapped square antiprismatic (CN=10)",
            "PCPA:11": "Pentagonal - face capped pentagonal antiprismatic (CN=11)",
            "H:11": "Hendecahedral (CN=11)",
            "SH:11": "Sphenoid hendecahedral (CN=11)",
            "CO:11": "Cs - octahedral (CN=11)",
            "DI:11": "Diminished icosahedral (CN=12)",
            "I:12": "Icosahedral (CN=12)",
            "PBP: 12": "Pentagonal - face bicapped pentagonal prismatic (CN=12)",
            "TT:12": "Truncated tetrahedral (CN=12)",
            "C:12": "Cuboctahedral (CN=12)",
            "AC:12": "Anticuboctahedral (CN=12)",
            "SC:12": "Square cupola (CN=12)",
            "S:12": "Sphenomegacorona (CN=12)",
            "HP:12": "Hexagonal prismatic (CN=12)",
            "HA:12": "Hexagonal antiprismatic (CN=12)",
            "SH:13": "Square-face capped hexagonal prismatic (CN=13)",
            "H:5": "H:5",
            "1": "1-fold",
            "2": "2-fold",
            "3": "3-fold",
            "4": "4-fold",
            "5": "5-fold",
            "6": "6-fold",
            "7": "7-fold",
            "8": "8-fold",
            "9": "9-fold",
            "10": "10-fold",
            "11": "11-fold",
            "12": "12-fold",
            "13": "13-fold",
            "14": "14-fold",
            "15": "15-fold",
            "16": "16-fold",
            "17": "17-fold",
            "18": "18-fold",
            "19": "19-fold",
            "20": "20-fold",
            "21": "21-fold",
            "22": "22-fold",
            "23": "23-fold",
            "24": "24-fold",
            "25": "25-fold",
            "26": "26-fold",
            "27": "27-fold",
            "28": "28-fold",
            "29": "29-fold",
            "30": "30-fold",
        }
        for key, items in results_dict.items():
            assert Description._coordination_environment_to_text(key) == items

    def test_plot(self, describe_nacl, describe_nacl_spin, describe_nacl_all, describe_csh_all):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp0:
            filename_test = Path(tmp0) / "test.pdf"
            describe_nacl.plot_cohps(save=True, filename=filename_test, xlim=[-4, 4])
            assert Path(filename_test).exists()

        with tempfile.TemporaryDirectory() as tmp1:
            filename_test = Path(tmp1) / "test.pdf"
            describe_nacl_spin.plot_cohps(save=True, filename=filename_test, xlim=[-4, 4])
            assert Path(filename_test).exists()

        with tempfile.TemporaryDirectory() as tmp2:
            filename_test = Path(tmp2) / "test.pdf"
            describe_nacl_all.plot_cohps(save=True, filename=filename_test, xlim=[-4, 4])
            filename_test_1 = Path(tmp2) / "test-0.pdf"
            filename_test_2 = Path(tmp2) / "test-1.pdf"
            assert not Path(filename_test).exists()
            assert Path(filename_test_1).exists()
            assert Path(filename_test_2).exists()

        with tempfile.TemporaryDirectory() as tmp2:
            filename_test = str(Path(tmp2) / "test.pdf")
            describe_nacl_all.plot_cohps(save=True, filename=filename_test, xlim=[-4, 4])
            filename_test_1 = Path(tmp2) / "test-0.pdf"
            filename_test_2 = Path(tmp2) / "test-1.pdf"
            assert not Path(filename_test).exists()
            assert Path(filename_test_1).exists()
            assert Path(filename_test_2).exists()

        with tempfile.TemporaryDirectory() as tmp4:
            filename_test = str(Path(tmp4) / "test.pdf")
            describe_csh_all.plot_cohps(save=True, filename=filename_test, xlim=[-4, 4])
            assert Path(filename_test).exists()

    def test_write_description(
        self,
        describe_nacl,
        describe_nasi_madelung_all,
        describe_nasbf6,
        describe_nasbf6_anbd,
        describe_cdf_anbd,
        describe_nacl_nan,
    ):
        describe_nacl.write_description()
        describe_nasi_madelung_all.write_description()
        describe_nasbf6.write_description()
        describe_nasbf6_anbd.write_description()
        describe_cdf_anbd.write_description()
        describe_nacl_nan.write_description()

    def test_text(
        self,
        describe_cdf,
        describe_nacl,
        describe_nasbf6,
        describe_nasbf6_anbd,
        describe_nacl_nan,
        describe_cdf_anbd,
        describe_k3sb,
        describe_k3sb_all,
        describe_csh_all,
        describe_batio3_orb,
        describe_c_orb,
        describe_nasbf6_orb,
        describe_cdf_comp_range_coop,
        describe_nacl_comp_range_cobi,
    ):
        assert describe_cdf.text == [
            "The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.",
            "Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOHP: -0.62 eV, "
            "44.26 percent antibonding interaction below EFermi) bonds.",
        ]
        # assert describe_NaCl.text == [
        #     "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
        #     "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV, "
        #     "3.448 percent antibonding interaction below EFermi) bonds.",
        # ]
        assert describe_nacl.text == [
            "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV,"
            " 3.448 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_nasbf6.text == [
            "The compound NaSbF6 has 2 symmetry-independent cation(s) with relevant cation-anion interactions: "
            "Na1, Sb2.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-F (mean ICOHP: -0.61 eV, "
            "4.071 percent antibonding interaction below EFermi) bonds.",
            "Sb2 has an octahedral (CN=6) coordination environment. It has 6 Sb-F (mean ICOHP: -5.45 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_nasbf6_anbd.text == [
            "The compound NaSbF6 has 2 symmetry-independent cation(s) with relevant cation-anion "
            "interactions: Na1, Sb2.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-F (mean ICOHP: -0.61 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
            "Sb2 has an octahedral (CN=6) coordination environment. It has 6 Sb-F (mean ICOHP: -5.45 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_nacl_nan.text == [
            "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_cdf_anbd.text == [
            "The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.",
            "Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOHP: -0.62 eV, "
            "100.0 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_k3sb.text == [
            "The compound K3Sb has 2 symmetry-independent cation(s) with relevant cation-anion "
            "interactions: K1, K2.",
            "K1 has a 6-fold coordination environment. It has 6 K-Sb (mean ICOHP: -0.14 eV, "
            "2.299 percent antibonding interaction below EFermi) bonds.",
            "K2 has a 4-fold coordination environment. It has 4 K-Sb (mean ICOHP: -0.36 eV, "
            "4.969 percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_k3sb_all.text == [
            "The compound K3Sb has 3 symmetry-independent atoms(s) with relevant bonds: K1, K2, Sb4.",
            "K1 has a 14-fold coordination environment. It has 8 K-K (mean ICOHP: -0.37 eV, 17.544 percent "
            "antibonding interaction below EFermi), and 6 K-Sb (mean ICOHP: -0.14 eV, 2.299 percent "
            "antibonding interaction below EFermi) bonds.",
            "K2 has a 14-fold coordination environment. It has 10 K-K (mean ICOHP: -0.22 eV, 17.073 "
            "percent antibonding interaction below EFermi), and 4 K-Sb (mean ICOHP: -0.36 eV, 4.969 "
            "percent antibonding interaction below EFermi) bonds.",
            "Sb4 has a 14-fold coordination environment. It has 14 Sb-K (mean ICOHP: -0.27 eV, 3.731 "
            "percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_csh_all.text == [
            "The compound CsH has 1 symmetry-independent atoms(s) with relevant bonds: Cs1.",
            "Cs1 has a 18-fold coordination environment. It has 18 Cs-Cs (mean ICOHP: -0.49 eV, 18.741 "
            "percent antibonding interaction below EFermi) bonds.",
        ]
        assert describe_batio3_orb.text == [
            "The compound BaTiO3 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Ti2.",
            "Ti2 has an octahedral (CN=6) coordination environment. It has 6 Ti-O (mean ICOHP: -3.54 eV, "
            "1.092 percent antibonding interaction below EFermi) bonds.",
            "In the 6 Ti-O bonds, relative to the summed ICOHPs, the maximum bonding contribution is "
            "from the Ti(3dz2)-O(2pz) orbital, contributing 14.0 percent, whereas the maximum "
            "antibonding contribution is from the Ti(4s)-O(2s) orbital, contributing 20.0 percent.",
        ]
        assert describe_c_orb.text == [
            "The compound C has 1 symmetry-independent atoms(s) with relevant bonds: C1.",
            "C1 has a tetrahedral (CN=4) coordination environment. "
            "It has 4 C-C (mean ICOHP: -9.59 eV, 0.0 percent antibonding interaction below EFermi) bonds.",
            "In the 4 C-C bonds, relative to the summed ICOHPs, the maximum bonding contribution is "
            "from the C(2s)-C(2s) orbital, contributing 10.0 percent, whereas no significant "
            "antibonding contribution is found in this bond.",
        ]
        assert describe_nasbf6_orb.text == [
            "The compound NaSbF6 has 3 symmetry-independent atoms(s) with relevant bonds: Na1, Sb2, F3.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-F (mean ICOHP: -0.61 eV, "
            "4.071 percent antibonding interaction below EFermi) bonds.",
            "In the 6 Na-F bonds, relative to the summed ICOHPs, the maximum bonding contribution is "
            "from the Na(3s)-F(2s) orbital, contributing 67.0 percent, whereas the maximum antibonding "
            "contribution is from Na(2py)-F(2s), Na(2pz)-F(2s), and Na(2px)-F(2s) orbitals, "
            "contributing 13.0, 13.0, and 13.0 percent, respectively.",
            "Sb2 has an octahedral (CN=6) coordination environment. It has 6 Sb-F (mean ICOHP: -5.45 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
            "In the 6 Sb-F bonds, relative to the summed ICOHPs, the maximum bonding contribution "
            "is from Sb(5py)-F(2s), Sb(5pz)-F(2s), and Sb(5px)-F(2s) orbitals, contributing 14.0, "
            "14.0, and 14.0 percent, respectively, whereas no significant antibonding contribution is "
            "found in this bond.",
            "F3 has a linear (CN=2) coordination environment. It has 1 F-Na (mean ICOHP: -0.61 eV, "
            "4.545 percent antibonding interaction below EFermi), and 1 F-Sb (mean ICOHP: -5.45 eV, "
            "0.0 percent antibonding interaction below EFermi) bonds.",
            "In the 1 F-Na bond, relative to the summed ICOHPs, the maximum bonding contribution "
            "is from the F(2s)-Na(3s) orbital, contributing 68.0 percent, whereas the maximum "
            "antibonding contribution is from F(2s)-Na(2pz) and F(2pz)-Na(2pz) orbitals, "
            "contributing 36.0 and 36.0 percent, respectively. In the 1 F-Sb bond, relative "
            "to the summed ICOHPs, the maximum bonding contribution is from the "
            "F(2s)-Sb(5pz) orbital, contributing 41.0 percent, whereas no significant "
            "antibonding contribution is found in this bond.",
        ]
        assert describe_cdf_comp_range_coop.text == [
            "The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.",
            "Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOOP: 0.01, 40.984 percent "
            "antibonding interaction below EFermi) bonds.",
        ]
        assert describe_nacl_comp_range_cobi.text == [
            "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
            "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOBI: 0.08, 0.0 percent "
            "antibonding interaction below EFermi) bonds.",
        ]


class TestCalcQualityDescribe:
    def test_calc_quality_description_text(self):
        calc_quality_K3Sb = Analysis.get_lobster_calc_quality_summary(
            path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
            path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
            path_to_lobsterout=TestDir / "test_data/K3Sb/lobsterout.gz",
            path_to_lobsterin=TestDir / "test_data/K3Sb/lobsterin.gz",
            potcar_symbols=["K_sv", "Sb"],
            path_to_bandoverlaps=TestDir / "test_data/K3Sb/bandOverlaps.lobster.gz",
            dos_comparison=True,
            bva_comp=True,
            path_to_doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
            e_range=[-20, 0],
            path_to_vasprun=TestDir / "test_data/K3Sb/vasprun.xml.gz",
            n_bins=256,
        )

        calc_quality_CsH = Analysis.get_lobster_calc_quality_summary(
            path_to_poscar=TestDir / "test_data/CsH/POSCAR.gz",
            path_to_charge=TestDir / "test_data/CsH/CHARGE.lobster.gz",
            path_to_lobsterout=TestDir / "test_data/CsH/lobsterout.gz",
            path_to_lobsterin=TestDir / "test_data/CsH/lobsterin.gz",
            potcar_symbols=["Cs_sv", "H"],
            path_to_bandoverlaps=TestDir / "test_data/CsH/bandOverlaps.lobster.gz",
            dos_comparison=False,
            bva_comp=True,
        )

        calc_quality_k3sb_des = Description.get_calc_quality_description(calc_quality_K3Sb)
        assert calc_quality_k3sb_des == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 0.83 and 6.36 %, respectively.",
            "The bandOverlaps.lobster file is generated during the LOBSTER run. This indicates that the "
            "projected wave function is not completely orthonormalized; however, the maximal deviation values "
            "observed compared to the identity matrix is below the threshold of 0.1.",
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis.",
            "The Tanimoto index from DOS comparisons in the energy range between -20, 0 eV for s, p,"
            " summed orbitals are: 0.8532, 0.9481, 0.9275.",
        ]

        calc_quality_csh_des = Description.get_calc_quality_description(calc_quality_CsH)
        assert calc_quality_csh_des == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 3.01 and 13.73 %, respectively.",
            "The bandOverlaps.lobster file is generated during the LOBSTER run. This indicates that the projected "
            "wave function is not completely orthonormalized. "
            "The maximal deviation value from the identity matrix is 0.4285, and there are 0.1822 percent "
            "k-points above the deviation threshold of 0.1. Please check the results of other quality checks "
            "like dos comparisons, charges, charge spillings before using the results for further analysis.",
            "The atomic charge signs from Mulliken population analysis agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis agree with the bond valence analysis.",
        ]


class TestCalcQualityDescribeWarnings:
    def test_warnings(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("once")
            calc_quality_warnings = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=TestDir / "test_data/BaTe_low_quality/POSCAR.gz",
                path_to_charge=TestDir / "test_data/BaTe_low_quality/CHARGE.lobster.gz",
                path_to_lobsterout=TestDir / "test_data/BaTe_low_quality/lobsterout.gz",
                path_to_lobsterin=TestDir / "test_data/BaTe_low_quality/lobsterin.gz",
                potcar_symbols=["Ba_sv", "Te"],
                path_to_doscar=TestDir / "test_data/BaTe_low_quality/DOSCAR.lobster.gz",
                path_to_vasprun=TestDir / "test_data/BaTe_low_quality/vasprun.xml.gz",
                e_range=[-50, 60],
                dos_comparison=True,
                bva_comp=False,
            )
        assert ("Consider using DOSCAR.LSO.lobster" in str(w[0].message))
        assert ("Minimum energy range requested" in str(w[1].message))
        assert ("Maximum energy range requested" in str(w[2].message))
        assert ("Input DOS files have very few points" in str(w[3].message))

        calc_des = Description.get_calc_quality_description(calc_quality_warnings)

        assert calc_des == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 2.255 and 12.72 %, respectively.",
            "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
            "generated during the LOBSTER run.",
            "The Tanimoto index from DOS comparisons in the energy range between -5, 0 eV for s, p, summed orbitals "
            "are: 0.4057, 0.2831, 0.2762.",
        ]

        with warnings.catch_warnings(record=True) as w2:
            warnings.simplefilter("once")
            calc_quality_warnings2 = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=TestDir / "test_data/C/POSCAR.gz",
                path_to_charge=TestDir / "test_data/C/CHARGE.lobster.gz",
                path_to_lobsterout=TestDir / "test_data/C/lobsterout.gz",
                path_to_lobsterin=TestDir / "test_data/C/lobsterin.gz",
                potcar_symbols=["C"],
                bva_comp=True,
            )
        assert ("Oxidation states from BVA analyzer cannot" in str(w2[0].message))

        calc_des2 = Description.get_calc_quality_description(calc_quality_warnings2)

        assert calc_des2 == [
            "The LOBSTER calculation used minimal basis.",
            "The absolute and total charge spilling for the calculation is 0.98 and 8.93 %, respectively.",
            "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is "
            "generated during the LOBSTER run.",
            "Oxidation states from BVA analyzer cannot be determined. Thus BVA charge comparison is not conducted.",
        ]

        with warnings.catch_warnings(record=True) as w3:
            warnings.simplefilter("once")
            calc_quality_warnings3 = Analysis.get_lobster_calc_quality_summary(
                path_to_poscar=TestDir / "test_data/BeTe/POSCAR.gz",
                path_to_charge=TestDir / "test_data/BeTe/CHARGE.lobster.gz",
                path_to_lobsterout=TestDir / "test_data/BeTe/lobsterout.gz",
                path_to_lobsterin=TestDir / "test_data/BeTe/lobsterin.gz",
                potcar_symbols=["Be_sv", "Te"],
                bva_comp=True,
            )
        assert "Consider rerunning the calc with the minimum basis" in str(w3[0].message)

        calc_des3 = Description.get_calc_quality_description(calc_quality_warnings3)

        assert calc_des3 == [
            "Consider rerunning the calculation with the minimum basis as well. "
            "Choosing a larger basis set is only recommended if you see a significant improvement of "
            "the charge spilling.",
            "The absolute and total charge spilling for the calculation is 1.48 and 13.99 %, respectively.",
            "The projected wave function is completely orthonormalized as no bandOverlaps.lobster file is generated "
            "during the LOBSTER run.",
            "The atomic charge signs from Mulliken population analysis do not agree with the bond valence analysis.",
            "The atomic charge signs from Loewdin population analysis do not agree with the bond valence analysis.",
        ]
