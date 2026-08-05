from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.io.lobster import Charge, Doscar, Icohplist

from lobsterpy.coxx.analyze import Analysis
from lobsterpy.coxx.describe import Description
from lobsterpy.featurize.core import FeaturizeIcoxxlist

TestDir = Path(__file__).absolute().parent


# Fixtures for testing analyze module
@pytest.fixture
def analyse_nacl():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_comp_range():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_comp_range_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )


@pytest.fixture
def analyse_nacl_comp_range_orb_with_objs():
    structure_obj = Structure.from_file(TestDir / "test_data/NaCl_comp_range/CONTCAR.gz")
    charge_obj = Charge(filename=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz")
    completecoxx_obj = CompleteCohp.from_file(
        filename=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
        fmt="LOBSTER",
        structure_file=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
    )
    icoxxlist_obj = Icohplist(filename=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz")

    return Analysis(
        structure=structure_obj,
        completecoxx_obj=completecoxx_obj,
        icoxxlist_obj=icoxxlist_obj,
        charge_obj=charge_obj,
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )


@pytest.fixture
def analyse_nacl_comp_range_cobi():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture
def analyse_nacl_comp_range_cobi_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
        orbital_resolved=True,
    )


@pytest.fixture
def analyse_nacl_nan():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )


@pytest.fixture
def analyse_nacl_valences():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=None,
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        type_charge="Valences",
    )


@pytest.fixture
def analyse_nacl_madelung():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        madelung_path=TestDir / "test_data/NaCl/MadelungEnergies.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_bati03():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/BaTiO3/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTiO3/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTiO3/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_batao2n1():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/BaTaO2N1/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_bati03_differentcutoff():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/BaTiO3/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTiO3/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTiO3/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.001,
    )


@pytest.fixture
def analyse_nacl_distorted():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_distorted/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_distorted/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_distorted/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_distorted/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_spin():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_spin/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_spin/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_spin/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_spin/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_all():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_madelung_all():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        madelung_path=TestDir / "test_data/NaCl/MadelungEnergies.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nasi_madelung_all():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaSi/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaSi/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaSi/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaSi/CHARGE.lobster.gz",
        madelung_path=TestDir / "test_data/NaSi/MadelungEnergies.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_batao2n1_cutoff():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/BaTaO2N1/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.001,
    )


@pytest.fixture
def analyse_nasbf6():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaSbF6/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nasbf6_anbd():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaSbF6/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-5.5,
    )


@pytest.fixture
def analyse_cdf():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/CdF/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )


@pytest.fixture
def analyse_cdf_comp_range():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_cdf_comp_range_coop():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF_comp_range/COOPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_coops=True,
    )


@pytest.fixture
def analyse_k3sb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/K3Sb/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_k3sb_all():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/K3Sb/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_k3sb_all_objs():
    structure_obj = Structure.from_file(TestDir / "test_data/K3Sb/CONTCAR.gz")
    charge_obj = Charge(filename=TestDir / "test_data/K3Sb/CHARGE.lobster.gz")
    completecoxx_obj = CompleteCohp.from_file(
        filename=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        fmt="LOBSTER",
        structure_file=TestDir / "test_data/K3Sb/CONTCAR.gz",
    )
    icoxxlist_obj = Icohplist(filename=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz")
    return Analysis(
        structure=structure_obj,
        completecoxx_obj=completecoxx_obj,
        charge_obj=charge_obj,
        icoxxlist_obj=icoxxlist_obj,
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_k3sb_all_cobi():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/K3Sb/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/K3Sb/COBICAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/K3Sb/ICOBILIST.lobster.gz",
        charge_path=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture
def analyse_k3sb_all_coop_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/K3Sb/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/K3Sb/COOPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/K3Sb/ICOOPLIST.lobster.gz",
        charge_path=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        orbital_resolved=True,
        are_coops=True,
    )


# fixtures for describe module tests
@pytest.fixture
def describe_cdf_comp_range_coop(analyse_cdf_comp_range_coop):
    return Description(analyse_cdf_comp_range_coop)


@pytest.fixture
def describe_nacl(analyse_nacl):
    return Description(analyse_nacl)


@pytest.fixture
def describe_nacl_valences(analyse_nacl_valences):
    return Description(analyse_nacl_valences)


@pytest.fixture
def describe_nacl_distorted(analyse_nacl_distorted):
    return Description(analyse_nacl_distorted)


@pytest.fixture
def describe_nacl_all(analyse_nacl_all):
    return Description(analyse_nacl_all)


@pytest.fixture
def describe_nacl_madelung_all(analyse_nacl_madelung_all):
    return Description(analyse_nacl_madelung_all)


@pytest.fixture
def describe_nacl_nan(analyse_nacl_nan):
    return Description(analyse_nacl_nan)


@pytest.fixture
def describe_nacl_comp_range_cobi(analyse_nacl_comp_range_cobi):
    return Description(analyse_nacl_comp_range_cobi)


@pytest.fixture
def describe_nasi_madelung_all(analyse_nasi_madelung_all):
    return Description(analyse_nasi_madelung_all)


@pytest.fixture
def describe_nasbf6(analyse_nasbf6):
    return Description(analyse_nasbf6)


@pytest.fixture
def describe_nasbf6_anbd(analyse_nasbf6_anbd):
    return Description(analyse_nasbf6_anbd)


@pytest.fixture
def describe_batao2n1():
    return Description(analyse_batao2n1)


# tests for key error issues
@pytest.fixture
def describe_k3sb(analyse_k3sb):
    return Description(analyse_k3sb)


@pytest.fixture
def describe_k3sb_all(analyse_k3sb_all):
    return Description(analyse_k3sb_all)


@pytest.fixture
def describe_batio3():
    analyse_batio3 = Analysis.from_files(
        structure_path=TestDir / "test_data/BaTiO3/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTiO3/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTiO3/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )
    return Description(analyse_batio3)


@pytest.fixture
def describe_batio3_orb():
    analyse_bati03_orb = Analysis.from_files(
        structure_path=TestDir / "test_data/BaTiO3/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTiO3/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTiO3/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )
    return Description(analyse_bati03_orb)


@pytest.fixture
def describe_c_orb():
    analyse_c_orb = Analysis.from_files(
        structure_path=TestDir / "test_data/C/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/C/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/C/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/C/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        orbital_resolved=True,
        summed_spins=False,
    )
    return Description(analyse_c_orb)


@pytest.fixture(scope="class")
def describe_cdf():
    analyse_cdf = Analysis.from_files(
        structure_path=TestDir / "test_data/CdF/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )
    return Description(analyse_cdf)


# test for empty bond dict text generation
@pytest.fixture
def describe_cdf_anbd():
    analyse_cdf_anbd = Analysis.from_files(
        structure_path=TestDir / "test_data/CdF/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )
    return Description(analyse_cdf_anbd)


@pytest.fixture
def describe_csh_all():
    analyse_csh_all = Analysis.from_files(
        structure_path=TestDir / "test_data/CsH/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CsH/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CsH/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CsH/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )
    return Description(analyse_csh_all)


@pytest.fixture
def describe_nacl_spin():
    analyse_nacl_spin = Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_spin/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_spin/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_spin/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_spin/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
    )
    return Description(analyse_nacl_spin)


@pytest.fixture
def describe_nasbf6_orb():
    analyse_nasbf6_orb = Analysis.from_files(
        structure_path=TestDir / "test_data/NaSbF6/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
        orbital_resolved=True,
    )
    return Description(analyse_nasbf6_orb)


# Fixtures for plotting module tests


@pytest.fixture
def bwdf_nacl():
    return FeaturizeIcoxxlist(
        path_to_icoxxlist=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        path_to_structure=TestDir / "test_data/NaCl/CONTCAR.gz",
        bin_width=0.1,
        normalization="formula_units",
        are_coops=False,
        are_cobis=False,
    ).calc_bwdf()


@pytest.fixture
def bwdf_cdf_cobi():
    return FeaturizeIcoxxlist(
        path_to_icoxxlist=TestDir / "test_data/CdF_comp_range/ICOBILIST.lobster.gz",
        path_to_structure=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
        bin_width=0.1,
        normalization="none",
        are_coops=False,
        are_cobis=True,
    ).calc_bwdf()


@pytest.fixture
def bwdf_cdf_coop():
    return FeaturizeIcoxxlist(
        path_to_icoxxlist=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
        path_to_structure=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
        bin_width=0.1,
        normalization="counts",
        are_coops=True,
        are_cobis=False,
    ).calc_site_bwdf(site_index=1)


@pytest.fixture
def plot_analyse_nacl():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
    )


@pytest.fixture
def plot_analyse_cdf_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/CdF_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        summed_spins=False,
        orbital_resolved=True,
    )


@pytest.fixture
def plot_analyse_nacl_cobi():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture
def plot_analyse_nacl_cobi_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        charge_path=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
        noise_cutoff=0.001,
        orbital_resolved=True,
        are_cobis=True,
    )


@pytest.fixture
def plot_analyse_nasi():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/NaSi/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/NaSi/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/NaSi/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/NaSi/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=True,
    )


@pytest.fixture
def plot_analyse_batio3_orb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/BaTiO3/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/BaTiO3/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/BaTiO3/CHARGE.lobster.gz",
        which_bonds="all",
        summed_spins=False,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )


@pytest.fixture
def plot_analyse_k3sb():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/K3Sb/CONTCAR.gz",
        coxxcar_path=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
    )


@pytest.fixture
def lobsterpy_plot_data():
    plot_data_file_name = TestDir / "test_data/interactive_plotter_ref/mp-8818.json.gz"

    with gzip.open(plot_data_file_name, "rb") as f:
        data = json.loads(f.read().decode("utf-8"))

    lobsterpy_plot_data = {}
    for item in data:
        lobsterpy_plot_data.update(item)

    return lobsterpy_plot_data["all_bonds"]["lobsterpy_data"]["cohp_plot_data"]


@pytest.fixture
def icohplist_nacl():
    return Icohplist(filename=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz")


@pytest.fixture
def icohplist_aln_v51():
    return Icohplist(filename=TestDir / "test_data/AlN_v51/ICOHPLIST.lobster.gz")


@pytest.fixture
def icooplist_nacl():
    return Icohplist(filename=TestDir / "test_data/NaCl_comp_range/ICOOPLIST.lobster.gz")


@pytest.fixture
def icobilist_nacl():
    return Icohplist(filename=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz")


@pytest.fixture
def nacl_dos():
    return Doscar(
        doscar=TestDir / "test_data/NaCl_comp_range/DOSCAR.LSO.lobster.gz",
        structure_file=TestDir / "test_data/NaCl_comp_range/CONTCAR.gz",
    )


@pytest.fixture
def k3sb_dos():
    return Doscar(
        doscar=TestDir / "test_data/K3Sb/DOSCAR.LSO.lobster.gz",
        structure_file=TestDir / "test_data/K3Sb/CONTCAR.gz",
    )


@pytest.fixture
def analyse_aln_v51():
    return Analysis.from_files(
        structure_path=TestDir / "test_data/AlN_v51/POSCAR.lobster.vasp.gz",
        coxxcar_path=TestDir / "test_data/AlN_v51/COHPCAR.lobster.gz",
        icoxxlist_path=TestDir / "test_data/AlN_v51/ICOHPLIST.lobster.gz",
        charge_path=TestDir / "test_data/AlN_v51/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )
