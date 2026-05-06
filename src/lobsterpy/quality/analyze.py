"""Analyze and summarize LOBSTER calculation quality."""

from __future__ import annotations

import warnings
from pathlib import Path

from monty.os.path import zpath
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.dos import CompleteDos, LobsterCompleteDos
from pymatgen.io.lobster import (
    Bandoverlaps,
    Charge,
    Doscar,
    Lobsterin,
    Lobsterout,
)
from pymatgen.io.vasp import Vasprun

from lobsterpy.quality.describe import get_calc_quality_description


class LobsterCalcQuality:
    """
    Analyze and summarize LOBSTER calculation quality.

    :param structure: pymatgen Structure object
    :param lobsterin: Lobsterin object
    :param lobsterout: Lobsterout object
    :param potcar_symbols: List of POTCAR symbols used in the calculation
    :param charge: Charge object (optional)
    :param bandoverlaps: Bandoverlaps object (optional)
    :param lobster_dos: LobsterCompleteDos object (optional)
    :param vasp_dos: CompleteDos object from VASP (optional)
    """

    def __init__(
        self,
        *,
        structure: Structure,
        lobsterin: Lobsterin,
        lobsterout: Lobsterout,
        potcar_symbols: list[str],
        charge: Charge | None = None,
        bandoverlaps: Bandoverlaps | None = None,
        lobster_dos: LobsterCompleteDos | None = None,
        vasp_dos: CompleteDos | None = None,
    ):
        """
        Initialize LobsterCalcQuality instance.

        :param structure: pymatgen Structure object
        :param lobsterin: Lobsterin object
        :param lobsterout: Lobsterout object
        :param potcar_symbols: List of POTCAR symbols used in the calculation
        :param charge: Charge object (optional)
        :param bandoverlaps: Bandoverlaps object (optional)
        :param lobster_dos: LobsterCompleteDos object (optional)
        :param vasp_dos: CompleteDos object from VASP (optional)
        """
        self.structure = structure
        self.lobsterin = lobsterin
        self.lobsterout = lobsterout
        self.potcar_symbols = potcar_symbols
        self.charge = charge
        self.bandoverlaps = bandoverlaps
        self.lobster_dos = lobster_dos
        self.vasp_dos = vasp_dos

    @classmethod
    def from_files(
        cls,
        *,
        poscar: Path,
        lobsterin: Path,
        lobsterout: Path,
        potcar: Path | None = None,
        potcar_symbols: list[str] | None = None,
        charge: Path | None = None,
        bandoverlaps: Path | None = None,
        doscar: Path | None = None,
        vasprun: Path | None = None,
    ) -> LobsterCalcQuality:
        """
        Create LobsterCalcQuality instance from explicit file paths.

        :param poscar: Path to CONTCAR/POSCAR file
        :param lobsterin: Path to lobsterin file
        :param lobsterout: Path to lobsterout file
        :param potcar: Path to POTCAR file (optional if potcar_symbols or vasprun provided)
        :param potcar_symbols: List of POTCAR symbols (optional if potcar or vasprun provided)
        :param charge: Path to CHARGE.lobster file (optional)
        :param bandoverlaps: Path to bandOverlaps.lobster file (optional)
        :param doscar: Path to DOSCAR.lobster file (optional)
        :param vasprun: Path to vasprun.xml file (optional)
        """
        structure = Structure.from_file(poscar)
        lobsterin_obj = Lobsterin.from_file(lobsterin)
        lobsterout_obj = Lobsterout(lobsterout)

        # POTCAR handling (logic unchanged)
        if potcar:
            potcar_symbols = Lobsterin._get_potcar_symbols(POTCAR_input=potcar)
        elif potcar_symbols:
            pass
        elif vasprun:
            vr = Vasprun(vasprun, parse_potcar_file=False, parse_eigen=False)
            potcar_symbols = [p.split()[1] for p in vr.potcar_symbols]
        else:
            raise ValueError(
                "Provide potcar, potcar_symbols, or vasprun.xml to identify basis set used for projections"
            )

        return cls(
            structure=structure,
            lobsterin=lobsterin_obj,
            lobsterout=lobsterout_obj,
            potcar_symbols=potcar_symbols,
            charge=Charge(filename=charge) if charge else None,
            bandoverlaps=Bandoverlaps(filename=bandoverlaps) if bandoverlaps and Path(bandoverlaps).exists() else None,
            lobster_dos=Doscar(doscar=doscar, structure_file=None, structure=structure).completedos if doscar else None,
            vasp_dos=Vasprun(vasprun, parse_potcar_file=False, parse_eigen=False).complete_dos if vasprun else None,
        )

    @classmethod
    def from_directory(cls, path_to_lobster_calc: Path | str) -> LobsterCalcQuality:
        """
        Create LobsterCalcQuality instance from a directory containing LOBSTER calculation files.

        :param path_to_lobster_calc: Path to directory with LOBSTER calculation files
        """
        if not isinstance(path_to_lobster_calc, Path):
            path_to_lobster_calc = Path(path_to_lobster_calc)

        return cls.from_files(
            poscar=zpath(path_to_lobster_calc / "CONTCAR"),
            lobsterin=zpath(path_to_lobster_calc / "lobsterin"),
            lobsterout=zpath(path_to_lobster_calc / "lobsterout"),
            potcar=zpath(path_to_lobster_calc / "POTCAR")
            if Path(zpath(path_to_lobster_calc / "POTCAR")).exists()
            else None,
            charge=zpath(path_to_lobster_calc / "CHARGE.lobster")
            if Path(zpath(path_to_lobster_calc / "CHARGE.lobster")).exists()
            else None,
            bandoverlaps=zpath(path_to_lobster_calc / "bandOverlaps.lobster")
            if Path(zpath(path_to_lobster_calc / "bandOverlaps.lobster")).exists()
            else None,
            doscar=zpath(path_to_lobster_calc / "DOSCAR.LSO.lobster")
            if Path(zpath(path_to_lobster_calc / "DOSCAR.LSO.lobster")).exists()
            else (
                zpath(path_to_lobster_calc / "DOSCAR.lobster")
                if Path(zpath(path_to_lobster_calc / "DOSCAR.lobster")).exists()
                else None
            ),
            vasprun=zpath(path_to_lobster_calc / "vasprun.xml")
            if Path(zpath(path_to_lobster_calc / "vasprun.xml")).exists()
            else None,
        )

    def get_calculation_quality_summary(
        self,
        *,
        dos_comparison: bool = False,
        bva_comp: bool = False,
        e_range: list[int] = [-5, 0],
        n_bins: int | None = None,
    ) -> dict:
        """
        Get a summary of the LOBSTER calculation quality.

        :param dos_comparison: Whether to include DOS comparison analysis (requires DOS data)
        :param bva_comp: Whether to include bond valence analysis (requires CHARGE.l
        :param e_range: Energy range for DOS comparison (default: [-5, 0] eV)
        :param n_bins: Number of bins for DOS comparison (default: None, uses maximum
            available bins in the specified energy range)

        Returns:
            Dictionary summarizing the calculation quality
        """
        quality = {}
        quality.update(self._minimal_basis())
        quality.update(self._charge_spilling())
        quality.update(self._band_overlaps())

        if bva_comp:
            quality.update(self._bva_charge_comparison())

        if dos_comparison:
            quality.update(self._dos_comparison(e_range, n_bins))

        return quality

    def _minimal_basis(self) -> dict:
        ref_bases = Lobsterin.get_all_possible_basis_functions(
            structure=self.structure,
            potcar_symbols=self.potcar_symbols,
        )

        calc_basis = [" ".join(b.split()[1:]) for b in self.lobsterin["basisfunctions"]]

        minimal = calc_basis == list(ref_bases[0].values())

        if not minimal:
            warnings.warn(
                "Consider rerunning the calc with the minimum basis as well. "
                "Choosing a larger basis set is recommended if charge "
                "spilling significantly improves.",
                stacklevel=2,
            )

        return {"minimal_basis": minimal}

    def _charge_spilling(self) -> dict:
        lob = self.lobsterout
        return {
            "charge_spilling": {
                "abs_charge_spilling": round(sum(lob.charge_spilling) / 2 * 100, 4),
                "abs_total_spilling": round(sum(lob.total_spilling) / 2 * 100, 4),
            }
        }

    def _band_overlaps(self) -> dict:
        if not self.bandoverlaps:
            return {
                "band_overlaps_analysis": {
                    "file_exists": False,
                    "limit_maxDeviation": None,
                    "has_good_quality_maxDeviation": True,
                    "max_deviation": None,
                    "percent_kpoints_abv_limit": None,
                }
            }

        total_kpoints = None
        for line in self.lobsterout.warning_lines:
            if "k-points could not be orthonormalized" in line:
                total_kpoints = int(line.split()[2])

        dev_above = [d for d in self.bandoverlaps.max_deviation if d > 0.1]
        assert total_kpoints is not None  # mypy type checking fails if not asserted

        return {
            "band_overlaps_analysis": {
                "file_exists": True,
                "limit_maxDeviation": 0.1,
                "has_good_quality_maxDeviation": self.bandoverlaps.has_good_quality_maxDeviation(0.1),
                "max_deviation": round(max(self.bandoverlaps.max_deviation), 4),
                "percent_kpoints_abv_limit": round(len(dev_above) / total_kpoints * 100, 4),
            }
        }

    def _bva_charge_comparison(self) -> dict:
        if not self.charge:
            raise ValueError("BVA comparison requested but CHARGE.lobster not provided")

        try:
            bva = BVAnalyzer()
            bva_oxi = ["POS" if v >= 0 else "NEG" for v in bva.get_valences(self.structure)]
            mull = ["POS" if v >= 0 else "NEG" for v in self.charge.Mulliken]
            loew = ["POS" if v >= 0 else "NEG" for v in self.charge.Loewdin]

            return {
                "charge_comparisons": {
                    "bva_mulliken_agree": mull == bva_oxi,
                    "bva_loewdin_agree": loew == bva_oxi,
                }
            }

        except ValueError:
            warnings.warn(
                "Oxidation states from BVA analyzer cannot be determined. Thus BVA charge comparison will be skipped",
                stacklevel=2,
            )
            return {"charge_comparisons": {}}

    def _dos_comparison(self, e_range, n_bins) -> dict:
        if not self.lobster_dos or not self.vasp_dos:
            raise ValueError("DOS comparison requested but DOS data missing")

        dos_lob = self.lobster_dos
        dos_vasp = self.vasp_dos

        min_e = int(max(e_range[0], min(dos_vasp.energies), min(dos_lob.energies)))
        max_e = int(min(e_range[-1], max(dos_vasp.energies), max(dos_lob.energies)))

        minimum_n_bins = min(
            len(dos_vasp.energies[(dos_vasp.energies >= min_e) & (dos_vasp.energies <= max_e)]),
            len(dos_lob.energies[(dos_lob.energies >= min_e) & (dos_lob.energies <= max_e)]),
        )

        if n_bins is None or n_bins > minimum_n_bins:
            n_bins = minimum_n_bins

        fp_kwargs = dict(min_e=min_e, max_e=max_e, n_bins=n_bins, normalize=True)

        dos_comp = {}

        for orb in dos_lob.get_spd_dos():
            tani = dos_vasp.get_dos_fp_similarity(
                dos_lob.get_dos_fp(fp_type=orb.name, **fp_kwargs),
                dos_vasp.get_dos_fp(fp_type=orb.name, **fp_kwargs),
                metric="tanimoto",
            )
            dos_comp[f"tanimoto_orb_{orb.name}"] = round(tani, 4)

        tani_sum = dos_vasp.get_dos_fp_similarity(
            dos_lob.get_dos_fp(fp_type="summed_pdos", **fp_kwargs),
            dos_vasp.get_dos_fp(fp_type="summed_pdos", **fp_kwargs),
            metric="tanimoto",
        )

        dos_comp.update(
            {
                "tanimoto_summed": round(tani_sum, 4),
                "e_range": [min_e, max_e],
                "n_bins": n_bins,
            }
        )

        return {"dos_comparisons": dos_comp}

    @staticmethod
    def describe(quality_dict: dict) -> list[str]:
        """
        Generate a text description of the LOBSTER calculation quality.

        :param quality_dict: python dictionary obtained from LobsterCalcQuality.get_calculation_quality_summary method

        Returns:
            list of strings describing the calculation quality

        """
        return get_calc_quality_description(quality_dict)

    @staticmethod
    def print_description(text: list[str]):
        """
        Print the calculation quality description to the screen.

        :param text: list of strings obtained from LobsterCalcQuality.describe method
        """
        print(" ".join(text))
