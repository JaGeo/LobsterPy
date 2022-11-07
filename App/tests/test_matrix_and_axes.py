import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

import numpy as np

from helper import get_coord_transformation_matrix


"""
def test_transformation_matrix():
    a = 1.0
    b = 2.0
    c = 3.0
    alpha = np.deg2rad(269.0)
    beta = np.deg2rad(-89.0)
    gamma = np.deg2rad(179.0)

    true_matrix = np.zeros(shape=(3,3))
    true_matrix[0][0] = a
    true_matrix[0][1] = b * np.cos(gamma)
    true_matrix[1][1] = b * np.sin(gamma)

    test_matrix = get_coord_transformation_matrix(a, b, c, alpha, beta, gamma)
    tol = 0.01

    diff00 = 0
    diff01 = 0
    diff02 = 0
    diff11 = 0
    diff12 = 0
    diff22 = 0



def test_axes():
    pass
"""
