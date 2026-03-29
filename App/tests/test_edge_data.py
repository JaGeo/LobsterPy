import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_coord_transformation_matrix, get_primitive_cell, get_primitive_supercell
from helpers_for_tests import lobstergraph, completecohp



def test_edge_properties_cell():
    cell = get_primitive_cell(lobstergraph, completecohp)

    for test_edge, (_, _, true_edge) in zip(cell["edges"].values(), lobstergraph.sg.graph.edges.data()):
        if test_edge["bond_length"] != true_edge["bond_length"]:
            assert False
        elif test_edge["icobi"] != true_edge["ICOBI"]:
            assert False
        elif test_edge["icoop"] != true_edge["ICOOP"]:
            assert False
        elif test_edge["icohp"] != true_edge["ICOHP"]:
            assert False
        elif test_edge["icohp_bonding_perc"] != true_edge["ICOHP_bonding_perc"]:
            assert False

    assert True



def test_edge_properties_supercell():
    a = lobstergraph.sg.structure.lattice.a
    b = lobstergraph.sg.structure.lattice.b
    c = lobstergraph.sg.structure.lattice.c
    alpha = lobstergraph.sg.structure.lattice.alpha
    beta = lobstergraph.sg.structure.lattice.beta
    gamma = lobstergraph.sg.structure.lattice.gamma

    frac_to_cart_matrix = get_coord_transformation_matrix(
        a, b, c, alpha, beta, gamma
    )

    cell = get_primitive_cell(lobstergraph, completecohp)
    cell = get_primitive_supercell(lobstergraph, cell, frac_to_cart_matrix)

    for test_edge, (_, _, true_edge) in zip(cell["edges"].values(), lobstergraph.sg.graph.edges.data()):
        if test_edge["bond_length"] != true_edge["bond_length"]:
            assert False
        elif test_edge["icobi"] != true_edge["ICOBI"]:
            assert False
        elif test_edge["icoop"] != true_edge["ICOOP"]:
            assert False
        elif test_edge["icohp"] != true_edge["ICOHP"]:
            assert False
        elif test_edge["icohp_bonding_perc"] != true_edge["ICOHP_bonding_perc"]:
            assert False

    assert True



def test_cohp_plot_cell():
    cell = get_primitive_cell(lobstergraph, completecohp)
    
    for bond_label, edge in cell["edges"].items():
        cohp_data = completecohp.get_cohp_by_label(label=bond_label).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x_true = [-(spinup_cohps[j] + spindown_cohps[j]) for j, _ in enumerate(spinup_cohps)]
        y_true = [energies[j] - fermi_energy for j, _ in enumerate(energies)]
        
        x_test = edge["cohp_plot"][0]
        y_test = edge["cohp_plot"][1]

        for test_point, true_point in zip(x_test, x_true):
            if test_point != true_point:
                assert False
        for test_point, true_point in zip(y_test, y_true):
            if test_point != true_point:
                assert False

    assert True



def test_cohp_plot_supercell():
    a = lobstergraph.sg.structure.lattice.a
    b = lobstergraph.sg.structure.lattice.b
    c = lobstergraph.sg.structure.lattice.c
    alpha = lobstergraph.sg.structure.lattice.alpha
    beta = lobstergraph.sg.structure.lattice.beta
    gamma = lobstergraph.sg.structure.lattice.gamma

    frac_to_cart_matrix = get_coord_transformation_matrix(
        a, b, c, alpha, beta, gamma
    )

    cell = get_primitive_cell(lobstergraph, completecohp)
    cell = get_primitive_supercell(lobstergraph, cell, frac_to_cart_matrix)

    for bond_label, edge in cell["edges"].items():
        cohp_data = completecohp.get_cohp_by_label(label=bond_label).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x_true = [-(spinup_cohps[j] + spindown_cohps[j]) for j, _ in enumerate(spinup_cohps)]
        y_true = [energies[j] - fermi_energy for j, _ in enumerate(energies)]

        x_test = edge["cohp_plot"][0]
        y_test = edge["cohp_plot"][1]

        for test_point, true_point in zip(x_test, x_true):
            if test_point != true_point:
                assert False
        for test_point, true_point in zip(y_test, y_true):
            if test_point != true_point:
                assert False

    assert True
