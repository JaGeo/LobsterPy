from pathlib import Path

from pymatgen.core import Structure

from lobsterpy.featurize import get_structure_path

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


def test_get_structure_path():
    """
    Tests that in case of both LOBSTER and VASP structure files present, the VASP
    file is read.

    :param lobster_path: path to root LOBSTER calc directory
    :return: path to structure file
    """
    poscar_path = get_structure_path(lobster_path=TestDir / "test_data/test_structure_path_handling")
    elements = Structure.from_file(poscar_path).elements
    assert "Zn" not in [el.symbol for el in elements]
