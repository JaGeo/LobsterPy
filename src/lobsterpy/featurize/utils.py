# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This package provides the modules for featurzing Lobster data ready for ML studies."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple
from warnings import warn

try:
    from mendeleev import element
except ImportError:
    element = None

import numpy as np
from monty.dev import requires
from monty.os.path import zpath

POSCAR_WARNING = (
    "Falling back to POSCAR, translations between individual atoms may differ from LOBSTER outputs. "
    "Please note that translations in the LOBSTER outputs are consistent with CONTCAR "
    "(also with POSCAR.lobster.vasp or POSCAR.vasp : written by LOBSTER >=v5)."
)


class CoxxFingerprint(NamedTuple):
    """
    Represents a Coxx fingerprint.

    This named tuple is used to store information related to a Coxx fingerprint, which
    includes energies, Coxx values, fingerprint type, spin type, number of bins, and bin width.

    :param energies: The energy values associated with the Coxx fingerprint.
    :param coxx: The Coxx values corresponding to each energy.
    :param fp_type: The type of the Coxx fingerprint.
    :param spin_type: The spin type associated with the fingerprint.
    :param n_bins: The number of bins used in the Coxx fingerprint.
    :param bin_width: The width of each bin in the Coxx fingerprint.
    """

    energies: np.ndarray
    coxx: np.ndarray
    fp_type: str
    spin_type: str
    n_bins: int
    bin_width: float


def get_file_paths(
    path_to_lobster_calc: str | Path = "", requested_files: list[str] = [], use_lso_dos: bool = True
) -> dict:
    """
    Get file paths for LobsterPy featurizations, raise Exception if not all of requested paths exist.

    :param path_to_lobster_calc: path to root LOBSTER calc directory
    :param requested_files: files to return paths for.
    :param use_lso_dos: solely required for BatchDosFeaturizer.
        Will force featurizer to use DOSCAR.LSO.lobster instead of DOSCAR.lobster.

    :return: dict that assigns each item of requested_files its path

    """
    default_values = {
        "structure": "CONTCAR",
        "cohpcar": "COHPCAR.lobster",
        "icohplist": "ICOHPLIST.lobster",
        "cobicar": "COBICAR.lobster",
        "icobilist": "ICOBILIST.lobster",
        "coopcar": "COOPCAR.lobster",
        "icooplist": "ICOOPLIST.lobster",
        "charge": "CHARGE.lobster",
        "madelung": "MadelungEnergies.lobster",
        "doscar": ("DOSCAR.LSO.lobster" if use_lso_dos else "DOSCAR.lobster"),
        "lobsterin": "lobsterin",
        "lobsterout": "lobsterout",
        "bandoverlaps": "bandOverlaps.lobster",
        "grosspop": "GROSSPOP.lobster",
        "potcar": "POTCAR",
        "vasprun": "vasprun.xml",
        "incar": "INCAR",
    }

    lobster_path = Path(path_to_lobster_calc)
    file_paths = {}
    missing_files = []

    for file in requested_files:
        file_str = default_values.get(file)
        file_str = file_str if isinstance(file_str, str) else file
        if file == "structure":
            try:
                file_paths[file] = get_structure_path(lobster_path=lobster_path)
            except Exception:
                missing_files.append(default_values["structure"])
        else:
            file_path = lobster_path / file_str
            if file_path.exists():
                file_paths[file] = file_path
            else:
                gz_file_path = Path(zpath(str(file_path.as_posix())))
                if gz_file_path.exists():
                    file_paths[file] = gz_file_path
                else:
                    missing_files.append(default_values[file])

    if missing_files:
        raise Exception(f"Files {missing_files} not found in {lobster_path.name}.")

    return file_paths


def get_structure_path(lobster_path: Path) -> Path:
    """
    Search iteratively for (unzipped / zipped) structure file.

    CONTCAR is prioritized over POSCAR.lobster.

    :param lobster_path: path to root LOBSTER calc directory

    :return: path to structure file
    """
    for filename in ["CONTCAR", "POSCAR.lobster", "POSCAR.lobster.vasp", "POSCAR"]:
        poscar_path = lobster_path / filename
        if poscar_path.exists():
            if filename == "POSCAR":
                warn(POSCAR_WARNING)
            return poscar_path
        gz_file_path = Path(zpath(str(poscar_path.as_posix())))
        if gz_file_path.exists():
            if filename == "POSCAR":
                warn(POSCAR_WARNING)
            return gz_file_path

    raise Exception


@requires(
    element is not None,
    "get_reduced_mass requires mendeleev. Reinstall package with `pip install lobsterpy[featurizer]`.",
)
def get_reduced_mass(atom_pair: list[str]) -> float:
    """
    Compute reduced mass between a pair of atoms.

    :param atom_pair: list of atomic species symbols in string

    :return: reduced mass
    """
    atom1 = element(atom_pair[0])
    atom2 = element(atom_pair[1])
    return (atom1.atomic_weight * atom2.atomic_weight) / (atom1.atomic_weight + atom2.atomic_weight)


@requires(
    element is not None,
    "get_electronegativities requires mendeleev. Reinstall package with `pip install lobsterpy[featurizer]`.",
)
def get_electronegativities(atom_pair: list[str]) -> list[float]:
    """
    Get Allen electronegativities for a pair of atoms.

    :param atom_pair: list of atomic species symbols in string

    :return: list of Allen electronegativities
    """
    atom1 = element(atom_pair[0])
    atom2 = element(atom_pair[1])
    return [atom1.electronegativity_allen(), atom2.electronegativity_allen()]


def sort_dict_by_value(input_dict: dict[str, float]) -> dict:
    """
    Sort dictionary by values.

    :param input_dict: input dictionary

    :return: sorted dictionary
    """
    return dict(sorted(input_dict.items(), key=lambda item: item[1]))
