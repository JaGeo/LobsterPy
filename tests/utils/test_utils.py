import gzip
import shutil
import warnings
from pathlib import Path

from pymatgen.core import Structure

from lobsterpy.utils import get_file_paths, get_structure_path

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


def test_get_structure_path_warning(tmp_path):
    # test that the warning is raised when reading POSCAR
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once")
        warnings.filterwarnings("ignore", module="pymatgen")
        source_file = TestDir / "test_data/test_structure_path_handling/CONTCAR.gz"
        temp_poscar_path = tmp_path / "POSCAR"  # copy CONTCAR as POSCAR
        shutil.copy(source_file, temp_poscar_path)
        _ = get_structure_path(lobster_path=tmp_path)

        assert (
            str(w[0].message) == "Falling back to POSCAR, translations between individual "
            "atoms may differ from LOBSTER outputs. Please note that "
            "translations in the LOBSTER outputs are consistent with "
            "CONTCAR (also with POSCAR.lobster.vasp or POSCAR.vasp : "
            "written by LOBSTER >=v5)."
        )

    # test that the warning is raised when reading POSCAR.gz
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("once")
        warnings.filterwarnings("ignore", module="pymatgen")
        source_file = TestDir / "test_data/test_structure_path_handling/CONTCAR.gz"
        temp_poscar_path = tmp_path / "POSCAR.gz"  # copy CONTCAR as POSCAR
        shutil.copy(source_file, temp_poscar_path)
        _ = get_structure_path(lobster_path=tmp_path)

        assert (
            str(w[0].message) == "Falling back to POSCAR, translations between individual "
            "atoms may differ from LOBSTER outputs. Please note that "
            "translations in the LOBSTER outputs are consistent with "
            "CONTCAR (also with POSCAR.lobster.vasp or POSCAR.vasp : "
            "written by LOBSTER >=v5)."
        )


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
