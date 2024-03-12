# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This package provides the modules for featurzing Lobster data ready for ML studies."""
from __future__ import annotations

from pathlib import Path

from monty.os.path import zpath


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
        "structure": "POSCAR",
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
                gz_file_path = Path(zpath(file_path))
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

    POSCAR is prioritized over POSCAR.lobster.

    :param lobster_path: path to root LOBSTER calc directory

    :return: path to structure file
    """
    for filename in ["POSCAR", "POSCAR.lobster"]:
        poscar_path = lobster_path / filename
        if poscar_path.exists():
            return poscar_path
        gz_file_path = Path(zpath(poscar_path))
        if gz_file_path.exists():
            return gz_file_path

    raise Exception
