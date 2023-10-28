from __future__ import annotations

import pytest
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description

CurrentDir = Path(__file__).absolute().parent
TestDir = CurrentDir / "../"


# Fixtures for testing analyze module
@pytest.fixture
def analyse_nacl():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_comp_range():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_comp_range_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )


@pytest.fixture
def analyse_nacl_comp_range_cobi():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture
def analyse_nacl_comp_range_cobi_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaCl_comp_range/COBICAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaCl_comp_range/ICOBILIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaCl_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
        orbital_resolved=True,
    )


@pytest.fixture
def analyse_nacl_nan():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )


@pytest.fixture
def analyse_nacl_valences():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=None,
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_madelung():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        path_to_madelung=TestDir / "test_data/NaCl/MadelungEnergies.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_bati03():
    return Analysis(
        path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
        path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_batao2n1():
    return Analysis(
        path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_bati03_differentcutoff():
    return Analysis(
        path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
        path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.001,
    )


@pytest.fixture
def analyse_nacl_distorted():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_distorted/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl_distorted/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl_distorted/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl_distorted/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_spin():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_spin/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl_spin/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl_spin/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl_spin/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_all():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nacl_madelung_all():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaCl/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl/CHARGE.lobster",
        path_to_madelung=TestDir / "test_data/NaCl/MadelungEnergies.lobster",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nasi_madelung_all():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaSi/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaSi/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaSi/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaSi/CHARGE.lobster",
        path_to_madelung=TestDir / "test_data/NaSi/MadelungEnergies.lobster",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_batao2n1_cutoff():
    return Analysis(
        path_to_poscar=TestDir / "test_data/BaTaO2N1/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/BaTaO2N1/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/BaTaO2N1/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/BaTaO2N1/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.001,
    )


@pytest.fixture
def analyse_nasbf6():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaSbF6/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_nasbf6_anbd():
    return Analysis(
        path_to_poscar=TestDir / "test_data/NaSbF6/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-5.5,
    )


@pytest.fixture
def analyse_cdf():
    return Analysis(
        path_to_poscar=TestDir / "test_data/CdF/POSCAR",
        path_to_cohpcar=TestDir / "test_data/CdF/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/CdF/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )


@pytest.fixture
def analyse_cdf_comp_range():
    return Analysis(
        path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_cdf_comp_range_coop():
    return Analysis(
        path_to_poscar=TestDir / "test_data/CdF_comp_range/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/CdF_comp_range/COOPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/CdF_comp_range/ICOOPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/CdF_comp_range/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_coops=True,
    )


@pytest.fixture
def analyse_k3sb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_k3sb_all():
    return Analysis(
        path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/K3Sb/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/K3Sb/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )


@pytest.fixture
def analyse_k3sb_all_cobi():
    return Analysis(
        path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/K3Sb/COBICAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/K3Sb/ICOBILIST.lobster.gz",
        path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        noise_cutoff=0.001,
        are_cobis=True,
    )


@pytest.fixture
def analyse_k3sb_all_coop_orb():
    return Analysis(
        path_to_poscar=TestDir / "test_data/K3Sb/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/K3Sb/COOPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/K3Sb/ICOOPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/K3Sb/CHARGE.lobster.gz",
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
    analyse_batio3 = Analysis(
        path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
        path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )
    return Description(analyse_batio3)


@pytest.fixture
def describe_batio3_orb():
    analyse_bati03_orb_orb = Analysis(
        path_to_poscar=TestDir / "test_data/BaTiO3/POSCAR",
        path_to_cohpcar=TestDir / "test_data/BaTiO3/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/BaTiO3/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/BaTiO3/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        orbital_cutoff=0.10,
        orbital_resolved=True,
    )
    return Description(analyse_bati03_orb_orb)


@pytest.fixture
def describe_c_orb():
    analyse_c_orb = Analysis(
        path_to_poscar=TestDir / "test_data/C/POSCAR",
        path_to_cohpcar=TestDir / "test_data/C/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/C/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/C/CHARGE.lobster",
        which_bonds="all",
        cutoff_icohp=0.1,
        orbital_resolved=True,
        summed_spins=False,
    )
    return Description(analyse_c_orb)


@pytest.fixture(scope="class")
def describe_cdf():
    analyse_cdf = Analysis(
        path_to_poscar=TestDir / "test_data/CdF/POSCAR",
        path_to_cohpcar=TestDir / "test_data/CdF/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/CdF/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
    )
    return Description(analyse_cdf)


# test for empty bond dict text generation
@pytest.fixture
def describe_cdf_anbd():
    analyse_cdf_anbd = Analysis(
        path_to_poscar=TestDir / "test_data/CdF/POSCAR",
        path_to_cohpcar=TestDir / "test_data/CdF/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/CdF/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/CdF/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        start=-4.0,
    )
    return Description(analyse_cdf_anbd)


@pytest.fixture
def describe_csh_all():
    analyse_csh_all = Analysis(
        path_to_poscar=TestDir / "test_data/CsH/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/CsH/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/CsH/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/CsH/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
    )
    return Description(analyse_csh_all)


@pytest.fixture
def describe_nacl_spin():
    analyse_nacl_spin = Analysis(
        path_to_poscar=TestDir / "test_data/NaCl_spin/POSCAR",
        path_to_cohpcar=TestDir / "test_data/NaCl_spin/COHPCAR.lobster",
        path_to_icohplist=TestDir / "test_data/NaCl_spin/ICOHPLIST.lobster",
        path_to_charge=TestDir / "test_data/NaCl_spin/CHARGE.lobster",
        which_bonds="cation-anion",
        cutoff_icohp=0.1,
        summed_spins=False,
    )
    return Description(analyse_nacl_spin)


@pytest.fixture
def describe_nasbf6_orb():
    analyse_nasbf6_orb = Analysis(
        path_to_poscar=TestDir / "test_data/NaSbF6/POSCAR.gz",
        path_to_cohpcar=TestDir / "test_data/NaSbF6/COHPCAR.lobster.gz",
        path_to_icohplist=TestDir / "test_data/NaSbF6/ICOHPLIST.lobster.gz",
        path_to_charge=TestDir / "test_data/NaSbF6/CHARGE.lobster.gz",
        which_bonds="all",
        cutoff_icohp=0.1,
        summed_spins=False,
        orbital_resolved=True,
    )
    return Description(analyse_nasbf6_orb)


# @pytest.fixture(scope="class")
# def describe_fixtures(
#     describe_nacl,
#     describe_nacl_nan,
#     describe_cdf,
#     describe_c_orb,
#     describe_nacl_all,
#     describe_k3sb_all,
#     describe_batio3,
#     describe_csh_all,
#     describe_batio3_orb,
#     describe_nasbf6,
#     describe_nasbf6_anbd,
#     describe_cdf_anbd,
#     describe_k3sb,
#     describe_nasbf6_orb,
#     describe_nasi_madelung_all,
#     describe_batao2n1,
#     describe_nacl_madelung_all,
#     describe_nacl_spin,
#     describe_nacl_distorted,
#     describe_nacl_valences,
#     describe_cdf_comp_range_coop,
#     describe_nacl_comp_range_cobi,
# ):
#     return {
#         "describe_nacl": describe_nacl,
#         "describe_nacl_nan": describe_nacl_nan,
#         "describe_cdf": describe_cdf,
#         "describe_c_orb": describe_c_orb,
#         "describe_nacl_all": describe_nacl_all,
#         "describe_k3sb_all": describe_k3sb_all,
#         "describe_batio3": describe_batio3,
#         "describe_csh_all": describe_csh_all,
#         "describe_batio3_orb": describe_batio3_orb,
#         "describe_nasbf6": describe_nasbf6,
#         "describe_nasbf6_anbd": describe_nasbf6_anbd,
#         "describe_cdf_anbd": describe_cdf_anbd,
#         "describe_k3sb": describe_k3sb,
#         "describe_nasbf6_orb": describe_nasbf6_orb,
#         "describe_nasi_madelung_all": describe_nasi_madelung_all,
#         "describe_batao2n1": describe_batao2n1,
#         "describe_nacl_madelung_all": describe_nacl_madelung_all,
#         "describe_nacl_spin": describe_nacl_spin,
#         "describe_nacl_distorted": describe_nacl_distorted,
#         "describe_nacl_valences": describe_nacl_valences,
#         "describe_cdf_comp_range_coop": describe_cdf_comp_range_coop,
#         "describe_nacl_comp_range_cobi": describe_nacl_comp_range_cobi,
#     }
