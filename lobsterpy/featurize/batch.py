# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines wrapper classes to quickly obtain similarity matrix of input fingerprint objects
"""
from __future__ import annotations
import os
from typing import NamedTuple
import multiprocessing as mp
from typing import List
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm
from lobsterpy.featurize.core import (
    FeaturizeLobsterpy,
    FeaturizeCharges,
    FeaturizeCOXX,
)

warnings.filterwarnings("ignore")


class BatchSummaryFeaturizer:
    """
    Featurizer sets that generates summary features from lobster data.

    Args:
        path_to_lobster_calcs: path to root directory consisting of all lobster calc
        path_to_jsons: path to root directory consisting of all lobster lightweight jsons
        feature_type: set the feature type for moment features.
        Possible options are "bonding", "antibonding" or "overall"
        charge_type : set charge type used for computing ionicity. Possible options are
        "mulliken", "loewdin or "both"
        bonds: "all_bonds" or "cation_anion_bonds"
        include_cobi_data : bool stating to include COBICAR.lobster features
        include_coop_data: bool stating to include COOPCAR.lobster features
        e_range : range of energy relative to fermi for which moment features needs to be computed
        n_jobs : parallel processes to run

    Attributes:
        get_df: a pandas dataframe with summary features
    """

    def __init__(
        self,
        path_to_lobster_calcs: str,
        path_to_jsons: str | None = None,
        feature_type: str = "antibonding",
        charge_type: str = "both",
        bonds: str = "all_bonds",
        include_cobi_data: bool = False,
        include_coop_data: bool = False,
        only_smallest_basis: bool = True,
        e_range: List[float] = [-5.0, 0.0],
        n_jobs: int = 4,
    ):
        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.path_to_jsons = path_to_jsons
        self.feature_type = feature_type
        self.charge_type = charge_type
        self.bonds = bonds
        self.include_cobi_data = include_cobi_data
        self.include_coop_data = include_coop_data
        self.only_smallest_basis = only_smallest_basis
        self.e_range = e_range
        self.n_jobs = n_jobs

    def _featurizelobsterpy(self, file_name_or_path):
        if os.path.isfile(file_name_or_path):
            featurize_lobsterpy = FeaturizeLobsterpy(
                path_to_json=file_name_or_path,
                bonds=self.bonds,
            )

        else:
            featurize_lobsterpy = FeaturizeLobsterpy(
                path_to_lobster_calc=file_name_or_path,
                bonds=self.bonds,
            )

        df = featurize_lobsterpy.get_df()

        return df

    def _featurizecoxx(self, path_to_lobster_calc):
        dir_name = Path(path_to_lobster_calc)
        coxxcar_path = dir_name / "COHPCAR.lobster.gz"
        structure_path = dir_name / "POSCAR.gz"
        icoxxlist_path = dir_name / "ICOHPLIST.lobster.gz"

        if coxxcar_path.exists() and structure_path.exists() and icoxxlist_path.exists():

            coxx = FeaturizeCOXX(
                path_to_coxxcar=coxxcar_path,
                path_to_icoxxlist=icoxxlist_path,
                path_to_structure=structure_path,
                feature_type=self.feature_type,
                e_range=self.e_range,
            )

            df_cohp = coxx.get_summarized_coxx_df()
        else:
            raise Exception("COHPCAR.lobster.gz or POSCAR.gz or ICOHPLIST.lobster.gz file "
                            "not found in {}".format(dir_name.name))

        if self.include_cobi_data:
            coxxcar_path = dir_name / "COBICAR.lobster.gz"
            icoxxlist_path = dir_name / "ICOBILIST.lobster.gz"

            if coxxcar_path.exists() and icoxxlist_path.exists():

                coxx = FeaturizeCOXX(
                    path_to_coxxcar=coxxcar_path,
                    path_to_icoxxlist=icoxxlist_path,
                    path_to_structure=structure_path,
                    feature_type=self.feature_type,
                    e_range=self.e_range,
                    are_cobis=True,
                )

                df_cobi = coxx.get_summarized_coxx_df()

            else:
                raise Exception("COBICAR.lobster.gz or ICOBILIST.lobster.gz file "
                                "not found in {}".format(dir_name.name))

        if self.include_coop_data:
            coxxcar_path = dir_name / "COOPCAR.lobster.gz"
            icoxxlist_path = dir_name / "ICOOPLIST.lobster.gz"

            if coxxcar_path.exists() and icoxxlist_path.exists():

                coxx = FeaturizeCOXX(
                    path_to_coxxcar=coxxcar_path,
                    path_to_icoxxlist=icoxxlist_path,
                    path_to_structure=structure_path,
                    feature_type=self.feature_type,
                    e_range=self.e_range,
                    are_coops=True,
                )

                df_coop = coxx.get_summarized_coxx_df()

            else:
                raise Exception("COOPCAR.lobster.gz or ICOOPLIST.lobster.gz file "
                                "not found in {}".format(dir_name.name))

        if self.include_cobi_data and self.include_coop_data:
            df = pd.concat([df_cohp, df_cobi, df_coop], axis=1)
        elif self.include_cobi_data and not self.include_coop_data:
            df = pd.concat([df_cohp, df_cobi], axis=1)
        elif not self.include_cobi_data and self.include_coop_data:
            df = pd.concat([df_cohp, df_coop], axis=1)
        else:
            df = df_cohp

        return df

    def _featurizecharges(self, path_to_lobster_calc):
        dir_name = Path(path_to_lobster_calc)
        charge_path = dir_name / "CHARGE.lobster.gz"
        structure_path = dir_name / "POSCAR.gz"

        if charge_path.exists() and structure_path.exists():
            if self.charge_type == "mulliken":
                charge_mull = FeaturizeCharges(
                    path_to_charge=charge_path,
                    path_to_structure=structure_path,
                    charge_type="mulliken",
                )
                df = charge_mull.get_df()
            elif self.charge_type == "loewdin":
                charge_loew = FeaturizeCharges(
                    path_to_charge=charge_path,
                    path_to_structure=structure_path,
                    charge_type="loewdin",
                )
                df = charge_loew.get_df()
            elif self.charge_type == "both":
                charge_mull = FeaturizeCharges(
                    path_to_charge=charge_path,
                    path_to_structure=structure_path,
                    charge_type="mulliken",
                )
                df_mull = charge_mull.get_df()

                charge_loew = FeaturizeCharges(
                    path_to_charge=charge_path,
                    path_to_structure=structure_path,
                    charge_type="loewdin",
                )
                df_loew = charge_loew.get_df()

                df = pd.concat([df_mull, df_loew], axis=1)

            return df
        else:
            raise Exception("CHARGE.lobster.gz or POSCAR.gz not found in {}".format(dir_name.name))


    def get_df(self):
        """
        This function featurizes LobsterPy condensed bonding analysis data from
        lobster lightweight json.gz files

        Returns:
            Returns a pandas dataframe with lobster summar icohp statistics

        """
        if self.path_to_jsons and self.only_smallest_basis:
            file_name_or_path = [
                os.path.join(self.path_to_jsons, f)
                for f in os.listdir(self.path_to_jsons)
                if not f.startswith("t")
                and not f.startswith(".")
                and not os.path.isdir(f)
                and not len(f.split("_")) > 1
            ]

        elif self.path_to_jsons and not self.only_smallest_basis:
            file_name_or_path = [
                os.path.join(self.path_to_jsons, f)
                for f in os.listdir(self.path_to_jsons)
                if not f.startswith("t")
                and not f.startswith(".")
                and not os.path.isdir(f)
            ]
        elif (
            self.path_to_lobster_calcs
            and self.only_smallest_basis
            and not self.path_to_jsons
        ):
            file_name_or_path = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
                and not len(f.split("_")) > 1
            ]

        elif (
            self.path_to_lobster_calcs
            and not self.only_smallest_basis
            and not self.path_to_jsons
        ):
            file_name_or_path = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
            ]

        with mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool:
            results = tqdm(
                pool.imap(self._featurizelobsterpy, file_name_or_path, chunksize=1),
                total=len(file_name_or_path),
                desc="Generating LobsterPy summary stats",
            )
            row = []
            for result in results:
                row.append(result)

        df_lobsterpy = pd.concat(row)

        if self.only_smallest_basis:
            paths = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
                and not len(f.split("_")) > 1
            ]
        else:
            paths = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
            ]

        with mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool:
            results = tqdm(
                pool.imap(self._featurizecoxx, paths, chunksize=1),
                total=len(paths),
                desc="Generating COHP/COOP/COBI summary stats",
            )
            row = []
            for result in results:
                row.append(result)

        df_coxx = pd.concat(row)

        with mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool:
            results = tqdm(
                pool.imap(self._featurizecharges, paths, chunksize=1),
                total=len(paths),
                desc="Generating charge based features",
            )
            row = []
            for result in results:
                row.append(result)

        df_charges = pd.concat(row)

        df = pd.concat([df_lobsterpy, df_coxx, df_charges], axis=1)

        return df


class BatchCoxxFingerprint:
    """
    Featurizer that generates Tanimoto index similarity matrix summary features from fingerprint objects.
    """

    def __init__(
        self,
        path_to_lobster_calcs: str,
        only_smallest_basis: bool = True,
        feature_type: str = "overall",
        label_list: List[str] | None = None,
        tanimoto: bool = True,
        normalize: bool = True,
        spin_type: str = "summed",
        n_bins: int = 56,
        e_range: List[float] = [-15.0, 0.0],
        n_jobs=4,
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.only_smallest_basis = only_smallest_basis
        self.feature_type = feature_type
        self.tanimoto = tanimoto
        self.normalize = normalize
        self.label_list = label_list
        self.spin_type = spin_type
        self.n_bins = n_bins
        self.e_range = e_range
        self.n_jobs = n_jobs
        self.are_cobis = are_cobis
        self.are_coops = are_coops

        self.fingerprint_df = self._get_fingerprints_df()

    def get_similarity_matrix_df(self):
        """
        This function will compute pairwise similarity index for each fingerprint object in input dataframe

        Returns:
             A Pandas dataframe
        """
        matrix = np.full(
            (self.fingerprint_df.shape[0], self.fingerprint_df.shape[0]), np.nan
        )
        for i, (row, col) in enumerate(self.fingerprint_df.iterrows()):
            for j, (row1, col1) in enumerate(self.fingerprint_df.iterrows()):
                #if i <= j:
                if self.tanimoto:
                    simi = self._get_fp_similarity(
                        col["COXX_FP"],
                        col1["COXX_FP"],
                        tanimoto=self.tanimoto,
                        normalize=False,
                    )
                else:
                    simi = self._get_fp_similarity(
                        col["COXX_FP"],
                        col1["COXX_FP"],
                        tanimoto=False,
                        normalize=True,
                    )
                matrix[i][j] = simi

        df = pd.DataFrame(
            matrix,
            index=list(self.fingerprint_df.index),
            columns=list(self.fingerprint_df.index),
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
        tanimoto: bool = True,
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
            BatchCoxxFingerprint._fp_to_dict(fp1) if not isinstance(fp1, dict) else fp1
        )

        fp2_dict = (
            BatchCoxxFingerprint._fp_to_dict(fp2) if not isinstance(fp2, dict) else fp2
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

    def _fingerprint_df(self, path_to_lobster_calc):
        dir_name = Path(path_to_lobster_calc)
        if self.are_cobis and not self.are_coops:
            if (dir_name / "COBICAR.lobster.gz").exists() and (dir_name / "ICOBILIST.lobster.gz").exists():
                coxxcar_path = dir_name / "COBICAR.lobster.gz"
                icoxxlist_path = dir_name / "ICOBILIST.lobster.gz"
            else:
                raise Exception("COBICAR.lobster.gz or ICOBILIST.lobster.gz file not found in "
                                "{}".format(dir_name.name))
        elif self.are_coops and not self.are_cobis:
            if (dir_name / "COOPCAR.lobster.gz").exists() and (dir_name / "ICOOPLIST.lobster.gz").exists():
                coxxcar_path = dir_name / "COOPCAR.lobster.gz"
                icoxxlist_path = dir_name / "ICOOPLIST.lobster.gz"
            else:
                raise Exception("COOPCAR.lobster.gz or ICOOPLIST.lobster.gz file not found in "
                                "{}".format(dir_name.name))
        else:
            if (dir_name / "COHPCAR.lobster.gz").exists() and (dir_name / "ICOHPLIST.lobster.gz").exists():
                coxxcar_path = dir_name / "COHPCAR.lobster.gz"
                icoxxlist_path = dir_name / "ICOHPLIST.lobster.gz"
            else:
                raise Exception("COHPCAR.lobster.gz or ICOHPLIST.lobster.gz file not found in "
                                "{}".format(dir_name.name))

        if (dir_name / "POSCAR.gz").exists():
            structure_path = dir_name / "POSCAR.gz"
        else:
            raise Exception("POSCAR.gz file not found in {}".format(dir_name.name))

        coxx = FeaturizeCOXX(
            path_to_coxxcar=coxxcar_path,
            path_to_icoxxlist=icoxxlist_path,
            path_to_structure=structure_path,
            feature_type=self.feature_type,
            e_range=self.e_range,
            are_coops=self.are_coops,
            are_cobis=self.are_cobis,
        )

        df_fp = coxx.get_coxx_fingerprint_df(
            spin_type=self.spin_type,
            n_bins=self.n_bins,
            normalize=self.normalize,
            label_list=self.label_list,
        )
        return df_fp

    def _get_fingerprints_df(self):
        if self.only_smallest_basis:
            paths = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
                and not len(f.split("_")) > 1
            ]
        else:
            paths = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
            ]

        with mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool:
            results = tqdm(
                pool.imap(self._fingerprint_df, paths, chunksize=1),
                total=len(paths),
                desc="Generating COHP/COBI/COOP fingerprints",
            )
            row = []
            for result in results:
                row.append(result)

        df = pd.concat(row)

        return df
