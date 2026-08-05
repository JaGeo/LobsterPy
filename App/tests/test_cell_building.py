import os
import sys

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

from helper import get_primitive_cell, get_primitive_supercell, get_coord_transformation_matrix
from helpers_for_tests import lobstergraph, completecohp
from itertools import product



def test_primitive_cell_building():
    cell = get_primitive_cell(lobstergraph, completecohp)

    true_vecs = [vec for vec in product([0, 1], repeat=3)]
    test_vecs = [cell["atoms"][bond_label]["frac_coord"] for bond_label in cell["atoms"].keys()]

    tol = 0.01

    for true_vec in true_vecs:
        flag = False
        for test_vec in test_vecs:
            diff_a = abs(true_vec[0] - test_vec[0])
            diff_b = abs(true_vec[1] - test_vec[1])
            diff_c = abs(true_vec[2] - test_vec[2])
            if (diff_a <= tol) and (diff_b <= tol) and (diff_c <= tol):
                flag = True
                break
        if not flag:
            assert False

    assert True



def test_primitive_cell_edges():
    cell = get_primitive_cell(lobstergraph, completecohp)

    tol = 0.01

    for edge in cell["edges"].values():
        for frac_coords in edge["frac_coords"]:
            start, end = frac_coords
            if not ((-tol <= start[0] <= 1+tol) and (-tol <= start[1] <= 1+tol) and (-tol <= start[2] <= 1+tol)):
                assert False
            if not ((-tol <= end[0] <= 1+tol) and (-tol <= end[1] <= 1+tol) and (-tol <= end[2] <= 1+tol)):
                assert False

    assert True



def test_primitive_supercell_building():
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

    data = list(lobstergraph.sg.graph.edges.data())

    xs = [vec["to_jimage"][0] for _, _, vec in data]
    a_true = 1 + abs(max(xs)) + abs(min(xs))
    x_max_test = 0
    x_min_test = 0

    ys = [vec["to_jimage"][1] for _, _, vec in data]
    b_true = 1 + abs(max(ys)) + abs(min(ys))
    y_max_test = 0
    y_min_test = 0

    zs = [vec["to_jimage"][2] for _, _, vec in data]
    c_true = 1 + abs(max(zs)) + abs(min(zs))
    z_max_test = 0
    z_min_test = 0

    test_vecs = [cell["atoms"][bond_label]["frac_coord"] for bond_label in cell["atoms"].keys()]

    tol = 0.01

    for test_vec in test_vecs:
        x_max_test = max(test_vec[0], x_max_test)
        x_min_test = min(test_vec[0], x_min_test)

        y_max_test = max(test_vec[1], y_max_test)
        y_min_test = min(test_vec[1], y_min_test)

        z_max_test = max(test_vec[2], z_max_test)
        z_min_test = min(test_vec[2], z_min_test)

    a_test = abs(x_max_test - x_min_test)
    b_test = abs(y_max_test - y_min_test)
    c_test = abs(z_max_test - z_min_test)

    if abs(a_true - a_test) <= tol and abs(b_true - b_test) <= tol and abs(c_true - c_test) <= tol:
        assert True
    else:
        assert False



def test_primitive_supercell_edges():
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

    data = list(lobstergraph.sg.graph.edges.data())

    xs = [vec["to_jimage"][0] for _, _, vec in data]
    x_max = 1 + max(xs)
    x_min = min(0, min(xs))

    ys = [vec["to_jimage"][1] for _, _, vec in data]
    y_max = 1 + max(ys)
    y_min = min(0, min(ys))

    zs = [vec["to_jimage"][2] for _, _, vec in data]
    z_max = 1 + max(zs)
    z_min = min(0, min(zs))

    tol = 0.01

    for edge in cell["edges"].values():
        for frac_coords in edge["frac_coords"]:
            start, end = frac_coords
            if not (
                (x_min - tol <= start[0] <= x_max + tol) and
                (y_min - tol <= start[1] <= y_max + tol) and
                (z_min - tol <= start[2] <= z_max + tol)
            ):
                assert False
            if not (
                (x_min - tol <= end[0] <= x_max + tol) and
                (y_min - tol <= end[1] <= y_max + tol) and
                (z_min - tol <= end[2] <= z_max + tol)
            ):
                assert False

    assert True
