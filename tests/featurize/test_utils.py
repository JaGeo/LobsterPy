from pathlib import Path

import pytest

from lobsterpy.featurize import (
    get_electronegativities,
    get_reduced_mass,
    sort_dict_by_value,
)

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


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
