# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This package provides the modules for featurzing Lobster data ready for ML studies."""

from __future__ import annotations

from typing import NamedTuple

try:
    from mendeleev import element
except ImportError:
    element = None

import numpy as np
from monty.dev import requires


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
