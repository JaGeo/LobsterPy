# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines wrapper classes to quickly obtain similarity matrix of input fingerprint objects
"""
from __future__ import annotations
from typing import NamedTuple
import numpy as np
import pandas as pd


class FingperintSimilarityMatrix:
    """
    Featurizer that generates Tanimoto index similarity matrix summary features from fingerprint objects.
    """

    def __init__(
        self,
        fingerpint_df: pd.DataFrame,
        tanimoto: bool = True,
        normalize: bool = False,
    ):
        self.fingerpint_df = fingerpint_df
        self.tanimoto = tanimoto
        self.normalize = normalize

    def get_similarity_matrix_df(self):
        """
        This function will compute pairwise similarity index for each fingerprint object in input dataframe

        Returns:
             A Pandas dataframe
        """
        matrix = np.full(
            (self.fingerpint_df.shape[0], self.fingerpint_df.shape[0]), np.nan
        )
        for i, (row, col) in enumerate(self.fingerpint_df.iterrows()):
            for j, (row1, col1) in enumerate(self.fingerpint_df.iterrows()):
                if i <= j:
                    simi = self._get_fp_similarity(
                        col["COXX_FP"],
                        col1["COXX_FP"],
                        tanimoto=self.tanimoto,
                        normalize=self.normalize,
                    )
                    matrix[i][j] = simi

        df = pd.DataFrame(
            matrix,
            index=list(self.fingerpint_df.index),
            columns=list(self.fingerpint_df.index),
        )

        return df

    @staticmethod
    def _fp_to_dict(fp) -> dict:
        """
        Converts a fingerprint into a dictionary

        Args:
            fp: The fingerprint to be converted into a dictionary

        Returns:
            dict: A dict of the fingerprint Keys=type, Values=np.ndarray(energies, cohp)
        """
        fp_dict = {}
        fp_dict[fp[2]] = np.array([fp[0], fp[1]], dtype="object").T

        return fp_dict

    @staticmethod
    def _get_fp_similarity(
        fp1,
        fp2,
        col: int = 1,
        pt: int | str = "All",
        normalize: bool = False,
        tanimoto: bool = False,
    ) -> float:
        """
        Calculates the similarity index (dot product) of two fingerprints

        Args:
            fp1 (NamedTuple): The 1st dos fingerprint object
            fp2 (NamedTuple): The 2nd dos fingerprint object
            col (int): The item in the fingerprints (0:energies,1: coxxs) to take the dot product of (default is 1)
            pt (int or str) : The index of the point that the dot product is to be taken (default is All)
            normalize (bool): If True normalize the scalar product to 1 (default is False)
            tanimoto (bool): If True will compute Tanimoto index (default is False)

        Raises:
            ValueError: If both tanimoto and normalize are set to True.

        Returns:
            Similarity index (float): The value of dot product
        """
        fp1_dict = (
            FingperintSimilarityMatrix._fp_to_dict(fp1)
            if not isinstance(fp1, dict)
            else fp1
        )

        fp2_dict = (
            FingperintSimilarityMatrix._fp_to_dict(fp2)
            if not isinstance(fp2, dict)
            else fp2
        )

        if pt == "All":
            vec1 = np.array([pt[col] for pt in fp1_dict.values()]).flatten()
            vec2 = np.array([pt[col] for pt in fp2_dict.values()]).flatten()
        else:
            vec1 = fp1_dict[fp1[2][pt]][col]
            vec2 = fp2_dict[fp2[2][pt]][col]

        if not normalize and tanimoto:
            rescale = (
                np.linalg.norm(vec1) ** 2
                + np.linalg.norm(vec2) ** 2
                - np.dot(vec1, vec2)
            )
            return np.dot(vec1, vec2) / rescale

        elif not tanimoto and normalize:
            rescale = np.linalg.norm(vec1) * np.linalg.norm(vec2)
            return np.dot(vec1, vec2) / rescale

        elif not tanimoto and not normalize:
            rescale = 1.0
            return np.dot(vec1, vec2) / rescale

        else:
            raise ValueError(
                "Cannot compute similarity index. Please set either normalize=True or tanimoto=True or both to False."
            )
