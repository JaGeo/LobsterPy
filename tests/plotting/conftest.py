from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from pymatgen.io.lobster import Doscar, Icohplist

from lobsterpy.cohp.analyze import Analysis

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


@pytest.fixture()
def plot_analyse_nacl():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
    )


@pytest.fixture()
def plot_analyse_cdf_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        summed_spins=False,
        orbital_resolved=True,
    )


@pytest.fixture()
def plot_analyse_nacl_cobi():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture()
def plot_analyse_nacl_cobi_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
        noise_cutoff=0.001,
        orbital_resolved=True,
        are_cobis=True,
    )


@pytest.fixture()
def plot_analyse_nasi():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaSi/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaSi/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaSi/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaSi/CHARGE.lobster",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=True,
    )


@pytest.fixture()
def plot_analyse_batio3_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
        path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
        which_bonds="all",
        summed_spins=False,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )


@pytest.fixture()
def plot_analyse_k3sb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
    )


@pytest.fixture()
def lobsterpy_plot_data():
    plot_data_file_name = TestDir / "test_data/interactive_plotter_ref/mp-8818.json.gz"

    with gzip.open(plot_data_file_name, "rb") as f:
        data = json.loads(f.read().decode("utf-8"))

    lobsterpy_plot_data = {}
    for item in data:
        lobsterpy_plot_data.update(item)

    return lobsterpy_plot_data["all_bonds"]["lobsterpy_data"]["cohp_plot_data"]


@pytest.fixture()
def icohplist_nacl():
    return Icohplist(
        filename=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz"
    )


@pytest.fixture()
def icooplist_nacl():
    return Icohplist(
        filename=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz"
    )


@pytest.fixture()
def icobilist_nacl():
    return Icohplist(
        filename=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz"
    )


@pytest.fixture()
def nacl_dos():
    return Doscar(
        doscar=TestDir / "test_data/NaCl_comp_range/DOSCAR.LSO.lobster.gz",
        structure_file=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
    )


@pytest.fixture()
def k3sb_dos():
    return Doscar(
        doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
        structure_file=TestDir / "test_data/K3Sb/POSCAR.gz",
    )
