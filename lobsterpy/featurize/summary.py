# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines wrapper class to obtain summary features ready to used for ML
"""

import logging
from typing import List
import pandas as pd
from lobsterpy.featurize.core import (
    FeaturizeLobsterpy,
    FeaturizeCharges,
    FeaturizeCOXX,
)


class SummarizedFeaturizer:
    """
    Featurizer sets that generates summary features from lobster data.

    Args:
        path_to_json: path to lobster lightweight json
        path_to_coxxcar: path to COXXCAR.lobster (e.g., "COXXCAR.lobster")
        path_to_icoxxlist : path to ICOXXLIST.lobster (e.g., "ICOXXLIST.lobster")
        path_to_structure : path to structure file (e.g., "POSCAR")
        feature_type: set the feature type for moment features and fingerprints.
        Possible options are "bonding", "antibonding" or "overall"
        path_to_charge : path to CHARGE.lobster (e.g., "CHARGE.lobster")
        charge_type : set charge type used for computing ionicity. Possible options are "Mulliken" or "Loewdin"
        bonds: "all_bonds" or "cation_anion_bonds"
        are_cobis : bool indicating if file contains COBI/ICOBI data
        are_coops : bool indicating if file contains COOP/ICOOP data
        e_range : range of energy relative to fermi for which moment features needs to be computed

    Attributes:
        get_summary_df: a pandas dataframe with summary features

    """

    def __init__(
        self,
        ids: str,
        path_to_json: str,
        path_to_coxxcar: str,
        path_to_icoxxlist: str,
        path_to_structure: str,
        feature_type: str,
        path_to_charge: str,
        charge_type: str,
        bonds: str = "all_bonds",
        e_range: List[float] = [-10.0, 0.0],
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
        self.ids = ids
        self.path_to_json = path_to_json
        self.bonds = bonds
        self.path_to_coxxcar = path_to_coxxcar
        self.path_to_icoxxlist = path_to_icoxxlist
        self.path_to_structure = path_to_structure
        self.feature_type = feature_type
        self.e_range = e_range
        self.are_cobis = are_cobis
        self.are_coops = are_coops
        self.path_to_charge = path_to_charge
        self.charge_type = charge_type

    def get_summary_df(self):
        """
        This function returns a pandas dataframe with all summary features column
        """
        featurize_lobsterpy = FeaturizeLobsterpy(
            path_to_json=self.path_to_json, bonds=self.bonds
        )

        df_lobsterpy = featurize_lobsterpy.get_df(ids=self.ids)

        coxx = FeaturizeCOXX(
            path_to_coxxcar=self.path_to_coxxcar,
            path_to_icoxxlist=self.path_to_icoxxlist,
            path_to_structure=self.path_to_structure,
            feature_type=self.feature_type,
            e_range=self.e_range,
        )

        df_coxx = coxx.get_summarized_coxx_df(ids=self.ids)

        feat_charge = FeaturizeCharges(
            path_to_charge=self.path_to_charge,
            path_to_structure=self.path_to_structure,
            charge_type=self.charge_type,
        )

        df_charge = feat_charge.get_df(ids=self.ids)

        df = pd.concat([df_lobsterpy, df_coxx, df_charge], axis=1)

        return df
