# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This package provides the module for featurzing LobsterPy lighweight jsons
"""

import gzip
import json
import numpy as np
import pandas as pd


def read_lobster_lightweight_json(filename: str):
    """
    This function reads loads the lightweight json.gz files and returns a python dictionary object
    with lobster summmarized bonding analysis data.
    Args:
        filename (str) : name of file
    Returns:
        Returns a dictionary with lobster summmarized bonding analysis data.
    """
    with gzip.open(filename, "rb") as f:
        data = json.loads(f.read().decode("utf-8"))

    lobster_data = {}
    for item in data:
        lobster_data.update(item)

    return lobster_data


def featurize_lobsterpy_icohp_data(filename, bonds: str = "all_bonds"):
    """
    This function reads loads the lightweight json.gz files and generates icohp stats data from lobsterpy
    bonding analysis
    Args:
        filename (str) : name of file
        bonds (str) : all_bonds | cation_anion_bonds
    Returns:
        Returns a pandas dataframe with lobsterpy icohp statistics .
    """
    # define a pandas dataframe
    df = pd.DataFrame(index=[filename.split(".")[0]])

    # read the lightweight lobster json files using read_lobster_lightweight_json function
    data = read_lobster_lightweight_json(filename)

    icohp_mean = []
    icohp_sum = []
    bond = []
    antibond = []
    # extract lobsterpy icohp related data for bond type specified
    for k, v in data[bonds]["lobsterpy_data"]["sites"].items():
        for k1, v1 in v["bonds"].items():
            icohp_mean.append(float(v1["ICOHP_mean"]))
            icohp_sum.append(float(v1["ICOHP_sum"]))
            bond.append(v1["bonding"]["perc"])
            antibond.append(v1["antibonding"]["perc"])

    # add stats data as columns to the dataframe
    df.loc[filename.split(".")[0], "Icohp_mean_avg"] = np.mean(icohp_mean)
    df.loc[filename.split(".")[0], "Icohp_mean_max"] = np.max(icohp_mean)
    df.loc[filename.split(".")[0], "Icohp_mean_min"] = np.min(icohp_mean)
    df.loc[filename.split(".")[0], "Icohp_mean_std"] = np.std(icohp_mean)

    df.loc[filename.split(".")[0], "Icohp_sum_avg"] = np.mean(icohp_sum)
    df.loc[filename.split(".")[0], "Icohp_sum_max"] = np.max(icohp_sum)
    df.loc[filename.split(".")[0], "Icohp_sum_min"] = np.min(icohp_sum)
    df.loc[filename.split(".")[0], "Icohp_sum_std"] = np.std(icohp_sum)

    df.loc[filename.split(".")[0], "bonding_perc_avg"] = np.mean(bond)
    df.loc[filename.split(".")[0], "bonding_perc_max"] = np.max(bond)
    df.loc[filename.split(".")[0], "bonding_perc_min"] = np.min(bond)
    df.loc[filename.split(".")[0], "bonding_perc_std"] = np.std(bond)

    df.loc[filename.split(".")[0], "antibonding_perc_avg"] = np.mean(antibond)
    df.loc[filename.split(".")[0], "antibonding_perc_min"] = np.min(antibond)
    df.loc[filename.split(".")[0], "antibonding_perc_max"] = np.max(antibond)
    df.loc[filename.split(".")[0], "antibonding_perc_std"] = np.std(antibond)

    # add madelung energies for the structure
    df.loc[filename.split(".")[0], "Madelung_Mull"] = data["madelung_energies"][
        "Mulliken"
    ]
    df.loc[filename.split(".")[0], "Madelung_Loew"] = data["madelung_energies"][
        "Loewdin"
    ]

    return df
