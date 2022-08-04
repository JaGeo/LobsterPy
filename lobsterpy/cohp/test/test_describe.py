import unittest
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../../"


# TODO: Add example without antibonding states


class TestDescribe(unittest.TestCase):
    def setUp(self):
        self.analyse_NaCl = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.describe_NaCl = Description(self.analyse_NaCl)

        self.analyse_NaCl_valences = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=None,
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.describe_NaCl_valences = Description(self.analyse_NaCl_valences)

        self.analyse_BaTiO3 = Analysis(
            path_to_poscar=TestDir / "TestData/BaTiO3/POSCAR",
            path_to_cohpcar=TestDir / "TestData/BaTiO3/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/BaTiO3/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/BaTiO3/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.describe_BaTiO3 = Description(self.analyse_BaTiO3)

        self.analyse_BaTaO2N1 = Analysis(
            path_to_poscar=TestDir / "TestData/BaTaO2N1/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/BaTaO2N1/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/BaTaO2N1/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/BaTaO2N1/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.describe_BaTaO2N1 = Description(self.analyse_BaTaO2N1)

        self.describe_BaTiO3 = Description(self.analyse_BaTiO3)

        self.analyse_CdF = Analysis(
            path_to_poscar=TestDir / "TestData/CdF/POSCAR",
            path_to_cohpcar=TestDir / "TestData/CdF/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/CdF/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/CdF/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.describe_CdF = Description(self.analyse_CdF)

        self.analyse_NaCl_distorted = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl_distorted/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl_distorted/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl_distorted/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl_distorted/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )

        self.describe_NaCl_distorted = Description(self.analyse_NaCl_distorted)

        self.analyse_NaCl_spin = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl_spin/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl_spin/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl_spin/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl_spin/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            summed_spins=False,
        )
        self.describe_NaCl_spin = Description(self.analyse_NaCl_spin)

        self.analyse_NaCl_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )

        self.describe_Nacl_all = Description(self.analyse_NaCl_all)

        self.analyse_NaCl_madelung_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            path_to_madelung=TestDir / "TestData/NaCl/MadelungEnergies.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )

        self.describe_Nacl_madelung_all = Description(self.analyse_NaCl_madelung_all)

        self.analyse_NaSi_madelung_all = Analysis(
            path_to_poscar=TestDir / "TestData/NaSi/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaSi/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaSi/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaSi/CHARGE.lobster",
            path_to_madelung=TestDir / "TestData/NaSi/MadelungEnergies.lobster",
            whichbonds="all",
            cutoff_icohp=0.1,
        )

        self.describe_NaSi_madelung_all = Description(self.analyse_NaSi_madelung_all)

        self.analyse_NaSbF6 = Analysis(
            path_to_poscar=TestDir / "TestData/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/NaSbF6/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
        )
        self.describe_NaSbF6 = Description(self.analyse_NaSbF6)

        self.analyse_NaSbF6_anbd = Analysis(
            path_to_poscar=TestDir / "TestData/NaSbF6/POSCAR.gz",
            path_to_cohpcar=TestDir / "TestData/NaSbF6/COHPCAR.lobster.gz",
            path_to_icohplist=TestDir / "TestData/NaSbF6/ICOHPLIST.lobster.gz",
            path_to_charge=TestDir / "TestData/NaSbF6/CHARGE.lobster.gz",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-5.5,
        )

        self.describe_NaSbF6_anbd = Description(self.analyse_NaSbF6_anbd)

        self.analyse_NaCl_nan = Analysis(
            path_to_poscar=TestDir / "TestData/NaCl/POSCAR",
            path_to_cohpcar=TestDir / "TestData/NaCl/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/NaCl/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/NaCl/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
        )

        self.describe_NaCl_nan = Description(self.analyse_NaCl_nan)

        self.analyse_CdF_anbd = Analysis(
            path_to_poscar=TestDir / "TestData/CdF/POSCAR",
            path_to_cohpcar=TestDir / "TestData/CdF/COHPCAR.lobster",
            path_to_icohplist=TestDir / "TestData/CdF/ICOHPLIST.lobster",
            path_to_charge=TestDir / "TestData/CdF/CHARGE.lobster",
            whichbonds="cation-anion",
            cutoff_icohp=0.1,
            start=-4.0,
        )

        self.describe_CdF_anbd = Description(self.analyse_CdF_anbd)

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
            "SA:8": "sqaure antiprismatic (CN=8)",
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
            self.assertEqual(Description._coordination_environment_to_text(key), items)

    def test_plot(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp0:
            filename_test = Path(tmp0) / "test.pdf"
            self.describe_NaCl.plot_cohps(
                save=True, filename=filename_test, xlim=[-4, 4]
            )
            self.assertTrue(Path(filename_test).exists())

        with tempfile.TemporaryDirectory() as tmp1:
            filename_test = Path(tmp1) / "test.pdf"
            self.describe_NaCl_spin.plot_cohps(
                save=True, filename=filename_test, xlim=[-4, 4]
            )
            self.assertTrue(Path(filename_test).exists())

        with tempfile.TemporaryDirectory() as tmp2:
            filename_test = Path(tmp2) / "test.pdf"
            self.describe_Nacl_all.plot_cohps(
                save=True, filename=filename_test, xlim=[-4, 4]
            )
            filename_test_1 = Path(tmp2) / "test-0.pdf"
            filename_test_2 = Path(tmp2) / "test-1.pdf"
            self.assertFalse(Path(filename_test).exists())
            self.assertTrue(Path(filename_test_1).exists())
            self.assertTrue(Path(filename_test_2).exists())

        with tempfile.TemporaryDirectory() as tmp2:
            filename_test = str(Path(tmp2) / "test.pdf")
            self.describe_Nacl_all.plot_cohps(
                save=True, filename=filename_test, xlim=[-4, 4]
            )
            filename_test_1 = Path(tmp2) / "test-0.pdf"
            filename_test_2 = Path(tmp2) / "test-1.pdf"
            self.assertFalse(Path(filename_test).exists())
            self.assertTrue(Path(filename_test_1).exists())
            self.assertTrue(Path(filename_test_2).exists())

    def test_write_descritoin(self):
        self.describe_NaCl.write_description()
        self.describe_NaSi_madelung_all.write_description()
        self.describe_NaSbF6.write_description()
        self.describe_NaSbF6_anbd.write_description()
        self.describe_CdF_anbd.write_description()
        self.describe_NaCl_nan.write_description()

    def test_text(self):
        self.assertEqual(
            self.describe_CdF.text,
            [
                "The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.",
                "Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOHP: -0.62 eV, 44.26 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_NaCl.text,
            [
                "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
                "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV, 3.448 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_NaCl.text,
            [
                "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
                "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV, 3.448 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_NaSbF6.text,
            [
                "The compound NaSbF6 has 2 symmetry-independent cation(s) with relevant cation-anion interactions: Na1, Sb2.",
                "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-F (mean ICOHP: -0.61 eV, 4.071 percent antibonding interaction below EFermi) bonds.",
                "Sb2 has an octahedral (CN=6) coordination environment. It has 6 Sb-F (mean ICOHP: -5.45 eV, 0.0 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_NaSbF6_anbd.text,
            [
                "The compound NaSbF6 has 2 symmetry-independent cation(s) with relevant cation-anion interactions: Na1, Sb2.",
                "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-F (mean ICOHP: -0.61 eV, 0.0 percent antibonding interaction below EFermi) bonds.",
                "Sb2 has an octahedral (CN=6) coordination environment. It has 6 Sb-F (mean ICOHP: -5.45 eV, 0.0 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_NaCl_nan.text,
            [
                "The compound NaCl has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Na1.",
                "Na1 has an octahedral (CN=6) coordination environment. It has 6 Na-Cl (mean ICOHP: -0.57 eV, 0.0 percent antibonding interaction below EFermi) bonds.",
            ],
        )
        self.assertEqual(
            self.describe_CdF_anbd.text,
            [
                "The compound CdF2 has 1 symmetry-independent cation(s) with relevant cation-anion interactions: Cd1.",
                "Cd1 has a cubic (CN=8) coordination environment. It has 8 Cd-F (mean ICOHP: -0.62 eV, 100.0 percent antibonding interaction below EFermi) bonds.",
            ],
        )


if __name__ == "__main__":
    unittest.main()
