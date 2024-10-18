import gzip
import shutil
from pathlib import Path

from pymatgen.core import Structure

from lobsterpy.featurize import get_file_paths, get_structure_path

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


def test_get_structure_path(tmp_path):
    """
    Tests that POSCAR path is returned by get_structure_path function.

    Tests that in case of both LOBSTER and VASP structure files present, the VASP
    file is read.
    """
    with (
        gzip.open(TestDir / "test_data/test_structure_path_handling/POSCAR.gz", "rb") as zipped_poscar,
        open(tmp_path / "POSCAR", "wb") as unzipped_poscar,
    ):
        shutil.copyfileobj(zipped_poscar, unzipped_poscar)

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
