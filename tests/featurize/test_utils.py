import gzip
import shutil
from pathlib import Path

import pytest
from pymatgen.core import Structure

from lobsterpy.featurize import (
    get_electronegativities,
    get_file_paths,
    get_reduced_mass,
    get_structure_path,
    sort_dict_by_value,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


def test_get_structure_path(tmp_path):
    """
    Tests that CONTCAR path is returned by get_structure_path function.

    Tests that in case of both LOBSTER and VASP structure files present, the VASP
    file is read.
    """
    with (
        gzip.open(TestDir / "test_data/test_structure_path_handling/CONTCAR.gz", "rb") as zipped_poscar,
        open(tmp_path / "CONTCAR", "wb") as unzipped_contcar,
    ):
        shutil.copyfileobj(zipped_poscar, unzipped_contcar)

    poscar_path_unzipped = get_structure_path(lobster_path=tmp_path)
    assert isinstance(poscar_path_unzipped, Path)

    poscar_path_both = get_structure_path(lobster_path=TestDir / "test_data/test_structure_path_handling")
    assert isinstance(poscar_path_both, Path)

    elements = Structure.from_file(poscar_path_both).elements
    assert "Zn" not in [el.symbol for el in elements]


def test_get_file_paths(tmp_path):
    """
    Tests that dict of str: Path is returned by get_file_paths().
    """
    file_paths_zipped = get_file_paths(
        path_to_lobster_calc=TestDir / "test_data/BaTaO2N1",
        requested_files=["structure", "cohpcar", "charge", "icohplist"],
    )
    for key, value in file_paths_zipped.items():
        assert isinstance(key, str)
        assert isinstance(value, Path)

    for file in ["COHPCAR.lobster", "ICOHPLIST.lobster"]:
        with (
            gzip.open(TestDir / f"test_data/BaTaO2N1/{file}.gz", "rb") as zipped_file,
            open(tmp_path / file, "wb") as unzipped_file,
        ):
            shutil.copyfileobj(zipped_file, unzipped_file)

    file_paths_unzipped = get_file_paths(path_to_lobster_calc=tmp_path, requested_files=["cohpcar", "icohplist"])
    for key, value in file_paths_unzipped.items():
        assert isinstance(key, str)
        assert isinstance(value, Path)


def test_get_reduced_mass():
    """
    Tests that reduced mass is computed correctly.
    """
    assert get_reduced_mass(["H", "Pt"]) == pytest.approx(1.002818, abs=1e-05)
    assert get_reduced_mass(["Na", "Cl"]) == pytest.approx(13.945765, abs=1e-05)


def test_get_electronegativities():
    """
    Tests that electronegativities are computed correctly.
    """
    assert get_electronegativities(["H", "Pt"]) == pytest.approx([13.61, 10.16], abs=1e-05)
    assert get_electronegativities(["Na", "Cl"]) == pytest.approx([5.14, 16.97], abs=1e-05)


def test_sort_dict_by_value():
    """
    Tests that dictionary is sorted by values.
    """
    assert sort_dict_by_value(input_dict={"Na-Cl": -1.6, "H-Pt": -3.4}) == {"H-Pt": -3.4, "Na-Cl": -1.6}
    assert sort_dict_by_value(input_dict={"a": 1, "b": 0, "c": 2}) == {"b": 0, "a": 1, "c": 2}
