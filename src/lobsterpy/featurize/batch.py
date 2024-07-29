# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module defines wrapper classes to quickly obtain similarity matrix of input fingerprint objects."""

from __future__ import annotations

import multiprocessing as mp
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm

from lobsterpy.featurize.core import (
    CoxxFingerprint,
    FeaturizeCharges,
    FeaturizeCOXX,
    FeaturizeDoscar,
    FeaturizeLobsterpy,
)
from lobsterpy.structuregraph.graph import LobsterGraph

from . import get_file_paths

warnings.filterwarnings("ignore")


class BatchSummaryFeaturizer:
    """
    Batch Featurizer sets that generates summary features from lobster data.

    :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
    :param path_to_jsons: path to root directory consisting of all lobster lightweight jsons
    :param feature_type: set the feature type for moment features.
        Possible options are `bonding`, `antibonding` or `overall`
    :param charge_type: set charge type used for computing ionicity. Possible options are
        `mulliken`, `loewdin` or `both`.
    :param bonds: `all_bonds` or `cation_anion_bonds`
    :param orbital_resolved: bool indicating whether LobsterPy analysis is performed orbital wise
    :param include_cobi_data: bool stating to include COBICAR.lobster features
    :param include_coop_data: bool stating to include COOPCAR.lobster features
    :param e_range: range of energy relative to fermi for which moment features needs to be computed
    :param n_jobs: parallel processes to run
    """

    def __init__(
        self,
        path_to_lobster_calcs: str | Path,
        path_to_jsons: str | Path | None = None,
        feature_type: str = "antibonding",
        charge_type: str = "both",
        bonds: str = "all",
        orbital_resolved: bool = False,
        include_cobi_data: bool = False,
        include_coop_data: bool = False,
        e_range: list[float] = [-5.0, 0.0],
        n_jobs: int = 4,
    ):
        """
        Featurize lobster data via multiprocessing for large number of compounds.

        :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
        :param path_to_jsons: path to root directory consisting of all lobster lightweight jsons
        :param feature_type: set the feature type for moment features.
            Possible options are `bonding`, `antibonding` or `overall`
        :param charge_type: set charge type used for computing ionicity. Possible options are
            `mulliken`, `loewdin` or `both`.
        :param bonds: `all` or `cation-anion` bonds
        :param orbital_resolved: bool indicating whether LobsterPy analysis is performed orbital wise
        :param include_cobi_data: bool stating to include COBICAR.lobster features
        :param include_coop_data: bool stating to include COOPCAR.lobster features
        :param e_range: range of energy relative to fermi for which moment features needs to be computed
        :param n_jobs: parallel processes to run
        """
        # Check for valid parameters of string type
        allowed_str_inputs = {
            "charge_type": ["mulliken", "loewdin", "both"],
            "bonds": ["all", "cation-anion"],
            "feature_type": ["bonding", "antibonding", "overall"],
        }
        for param, param_string in zip([charge_type, bonds, feature_type], ["charge_type", "bonds", "feature_type"]):
            if param not in allowed_str_inputs[param_string]:
                raise ValueError(
                    f"Parameter {param_string} set to {param} but must be in "
                    f"{list(allowed_str_inputs[param_string])}."
                )

        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.path_to_jsons = path_to_jsons
        self.feature_type = feature_type
        self.charge_type = charge_type
        self.bonds = bonds
        self.orbital_resolved = orbital_resolved
        self.include_cobi_data = include_cobi_data
        self.include_coop_data = include_coop_data
        self.e_range = e_range
        self.n_jobs = n_jobs

    def _featurizelobsterpy(self, file_name_or_path: str | Path) -> pd.DataFrame:
        """
        Featurize Lobsterpy condensed bonding analysis data.

        if lightweight json file exists loads that or invokes LobsterPy Analysis class.

        :param file_name_or_path: path to the LOBSTER calc directory or
            lightweight condensed bonding analysis json file name.

        Returns:
            A pandas dataframe with ICOHP stats like mean, min, max of relevant bonds and
            madelung energies

        """
        if Path(file_name_or_path).is_file():
            featurize_lobsterpy = FeaturizeLobsterpy(
                path_to_json=file_name_or_path,
                bonds=self.bonds,
            )

        else:
            featurize_lobsterpy = FeaturizeLobsterpy(
                path_to_lobster_calc=file_name_or_path,
                bonds=self.bonds,
                orbital_resolved=self.orbital_resolved,
            )

        return featurize_lobsterpy.get_df()

    def _featurizecoxx(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Featurize COHP/COBI/COOPCAR data using FeaturizeCOXX.

        :param path_to_lobster_calc: path to root LOBSTER calc directory

        Returns:
            A pandas dataframe with COHP summary stats data mainly weighted ICOHP/ICOOP/ICOBI,
            Effective interaction number and moment features (center, width, skewness and kurtosis)

        """
        file_paths = get_file_paths(
            path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "cohpcar", "icohplist"]
        )
        structure_path = file_paths.get("structure")

        coxx = FeaturizeCOXX(
            path_to_coxxcar=str(file_paths.get("cohpcar")),
            path_to_icoxxlist=str(file_paths.get("icohplist")),
            path_to_structure=str(structure_path),
            feature_type=self.feature_type,
            e_range=self.e_range,
        )

        df = coxx.get_summarized_coxx_df()
        del coxx

        if self.include_cobi_data:
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc, requested_files=["cobicar", "icobilist"]
            )

            coxx = FeaturizeCOXX(
                path_to_coxxcar=str(file_paths.get("cobicar")),
                path_to_icoxxlist=str(file_paths.get("icobilist")),
                path_to_structure=str(structure_path),
                feature_type=self.feature_type,
                e_range=self.e_range,
                are_cobis=True,
            )

            df_cobi = coxx.get_summarized_coxx_df()
            df = pd.concat([df, df_cobi], axis=1)
            del coxx

        if self.include_coop_data:
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc, requested_files=["coopcar", "icooplist"]
            )

            coxx = FeaturizeCOXX(
                path_to_coxxcar=str(file_paths.get("coopcar")),
                path_to_icoxxlist=str(file_paths.get("icooplist")),
                path_to_structure=str(structure_path),
                feature_type=self.feature_type,
                e_range=self.e_range,
                are_coops=True,
            )

            df_coop = coxx.get_summarized_coxx_df()
            df = pd.concat([df, df_coop], axis=1)
            del coxx

        return df

    def _featurizecharges(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Featurize CHARGE.lobster.gz data that using FeaturizeCharges.

        :param path_to_lobster_calc: path to root LOBSTER calc directory

        Returns:
            A pandas dataframe with computed ionicity for the structure

        """
        file_paths = get_file_paths(path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "charge"])

        if self.charge_type == "mulliken":
            charge_mull = FeaturizeCharges(
                path_to_charge=str(file_paths.get("charge")),
                path_to_structure=str(file_paths.get("structure")),
                charge_type="mulliken",
            )
            df = charge_mull.get_df()
        elif self.charge_type == "loewdin":
            charge_loew = FeaturizeCharges(
                path_to_charge=str(file_paths.get("charge")),
                path_to_structure=str(file_paths.get("structure")),
                charge_type="loewdin",
            )
            df = charge_loew.get_df()
        else:
            charge_mull = FeaturizeCharges(
                path_to_charge=str(file_paths.get("charge")),
                path_to_structure=str(file_paths.get("structure")),
                charge_type="mulliken",
            )
            df_mull = charge_mull.get_df()

            charge_loew = FeaturizeCharges(
                path_to_charge=str(file_paths.get("charge")),
                path_to_structure=str(file_paths.get("structure")),
                charge_type="loewdin",
            )
            df_loew = charge_loew.get_df()

            df = pd.concat([df_mull, df_loew], axis=1)

        return df

    def get_df(self) -> pd.DataFrame:
        """
        Generate a pandas dataframe with summary features extracted from LOBSTER files.

        Uses multiprocessing to speed up the process.

        Returns:
            Returns a pandas dataframe

        """
        if self.path_to_jsons:
            file_name_or_path = [
                os.path.join(self.path_to_jsons, f)
                for f in os.listdir(self.path_to_jsons)
                if not f.startswith("t") and not f.startswith(".") and not os.path.isdir(f)
            ]

        elif self.path_to_lobster_calcs and not self.path_to_jsons:
            file_name_or_path = [
                os.path.join(self.path_to_lobster_calcs, f)
                for f in os.listdir(self.path_to_lobster_calcs)
                if not f.startswith("t")
                and not f.startswith(".")
                and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
            ]

        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(file_name_or_path), desc="Generating LobsterPy summary stats") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._featurizelobsterpy, file_name_or_path, chunksize=1)):
                pbar.update()
                row.append(result)

        df_lobsterpy = pd.concat(row)
        df_lobsterpy.sort_index(inplace=True)  # noqa: PD002

        paths = [
            os.path.join(self.path_to_lobster_calcs, f)
            for f in os.listdir(self.path_to_lobster_calcs)
            if not f.startswith("t")
            and not f.startswith(".")
            and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
        ]

        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(paths), desc="Generating COHP/COOP/COBI summary stats") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._featurizecoxx, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df_coxx = pd.concat(row)
        df_coxx.sort_index(inplace=True)  # noqa: PD002

        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(paths), desc="Generating charge based features") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._featurizecharges, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df_charges = pd.concat(row)
        df_charges.sort_index(inplace=True)  # noqa: PD002

        return pd.concat([df_lobsterpy, df_coxx, df_charges], axis=1)


class BatchCoxxFingerprint:
    """
    BatchFeaturizer to generate COHP/COOP/COBI fingerprints and Tanimoto index similarity matrix.

    :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
    :param feature_type: set the feature type for moment features.
        Possible options are `bonding`, `antibonding` or `overall`
    :param label_list: bond labels list for which fingerprints needs to be generated.
    :param tanimoto: bool to state to compute tanimoto index between fingerprint objects
    :param normalize: bool to state to normalize the fingerprint data
    :param spin_type: can be `summed` or `up` or `down`.
    :param n_bins: sets number for bins for fingerprint objects
    :param e_range: range of energy relative to fermi for which moment features needs to be computed
    :param n_jobs: number of parallel processes to run
    :param fingerprint_for: Possible options are `cohp` or `cobi` or `coop`.
        Based on this fingerprints will be computed for COHPCAR/COOBICAR/COOPCAR.lobster files

    Attributes:
        fingerprint_df: A pandas dataframe with fingerprint objects
    """

    def __init__(
        self,
        path_to_lobster_calcs: str | Path,
        feature_type: str = "overall",
        label_list: list[str] | None = None,
        tanimoto: bool = True,
        normalize: bool = True,
        spin_type: str = "summed",
        n_bins: int = 56,
        e_range: list[float] = [-15.0, 0.0],
        n_jobs=4,
        fingerprint_for: str = "cohp",
    ):
        """
        Generate COHP/COOP/COBI fingerprints and pair-wise Tanimoto index similarity matrix.

        :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
        :param feature_type: set the feature type for moment features.
            Possible options are `bonding`, `antibonding` or `overall`
        :param label_list: bond labels list for which fingerprints needs to be generated.
        :param tanimoto: bool to state to compute tanimoto index between fingerprint objects
        :param normalize: bool to state to normalize the fingerprint data
        :param spin_type: can be `summed` or `up` or `down`.
        :param n_bins: sets number for bins for fingerprint objects
        :param e_range: range of energy relative to fermi for which moment features needs to be computed
        :param n_jobs: number of parallel processes to run
        :param fingerprint_for: Possible options are `cohp` or `cobi` or `coop`.
            Based on this fingerprints will be computed for COHPCAR/COOBICAR/COOPCAR.lobster files
        """
        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.feature_type = feature_type
        self.tanimoto = tanimoto
        self.normalize = normalize
        self.label_list = label_list
        self.spin_type = spin_type
        self.n_bins = n_bins
        self.e_range = e_range
        self.n_jobs = n_jobs
        self.fingerprint_for = fingerprint_for

        self.fingerprint_df = self._get_fingerprints_df()

    def get_similarity_matrix_df(self) -> pd.DataFrame:
        """
        Compute pairwise similarity index for each fingerprint object in input dataframe.

        Returns:
             A Pandas dataframe
        """
        matrix = np.full((self.fingerprint_df.shape[0], self.fingerprint_df.shape[0]), np.nan)
        for i, (_, col) in enumerate(self.fingerprint_df.iterrows()):
            for j, (_, col1) in enumerate(self.fingerprint_df.iterrows()):
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

        return pd.DataFrame(
            matrix,
            index=list(self.fingerprint_df.index),
            columns=list(self.fingerprint_df.index),
        )

    @staticmethod
    def _fp_to_dict(fp: CoxxFingerprint) -> dict:
        """
        Convert a fingerprint obj into a dictionary.

        :param fp: The fingerprint to be converted into a dictionary

        Returns:
            dict: A dict of the fingerprint Keys=type, Values=np.ndarray(energies, cohp)
        """
        fp_dict = {}
        fp_dict[fp[2]] = np.array([fp[0], fp[1]], dtype="object").T

        return fp_dict

    @staticmethod
    def _get_fp_similarity(
        fp1: CoxxFingerprint,
        fp2: CoxxFingerprint,
        col: int = 1,
        pt: int | str = "All",
        normalize: bool = False,
        tanimoto: bool = True,
    ) -> float:
        """
        Calculate the similarity index (dot product) of two fingerprints.

        :param fp1 The 1st CoxxFingerprint object
        :param fp2: The 2nd CoxxFingerprint object
        :param col: The item in the fingerprints (0:energies,1: coxxs) to take the dot product of (default is 1)
        :param pt: The index of the point that the dot product is to be taken (default is All)
        :param normalize: If True normalize the scalar product to 1 (default is False)
        :param tanimoto: If True will compute Tanimoto index (default is False)

        Raises:
            ValueError: If both tanimoto and normalize are set to True.

        Returns:
            Similarity index (float): The value of dot product

        """
        fp1_dict = BatchCoxxFingerprint._fp_to_dict(fp1) if not isinstance(fp1, dict) else fp1

        fp2_dict = BatchCoxxFingerprint._fp_to_dict(fp2) if not isinstance(fp2, dict) else fp2

        if pt == "All":
            vec1 = np.array([pt[col] for pt in fp1_dict.values()]).flatten()
            vec2 = np.array([pt[col] for pt in fp2_dict.values()]).flatten()
        else:
            vec1 = fp1_dict[fp1[2][pt]][col]  # type: ignore
            vec2 = fp2_dict[fp2[2][pt]][col]  # type: ignore

        if not normalize and tanimoto:
            rescale = np.linalg.norm(vec1) ** 2 + np.linalg.norm(vec2) ** 2 - np.dot(vec1, vec2)

        elif not tanimoto and normalize:
            rescale = np.linalg.norm(vec1) * np.linalg.norm(vec2)

        elif not tanimoto and not normalize:
            rescale = 1.0

        else:
            raise ValueError(
                "Cannot compute similarity index. Please set either normalize=True or tanimoto=True or both to False."
            )
        return np.dot(vec1, vec2) / rescale

    def _fingerprint_df(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Get fingerprint object dataframe via  FeaturizeCOXX.get_coxx_fingerprint_df.

        Also helps to generate the data used for fingerprint generation.

        :param path_to_lobster_calc: path to root LOBSTER calculation directory.

        Returns:
            A pandas dataframe with COXX fingerprint object

        """
        if self.fingerprint_for.upper() == "COBI":
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "cobicar", "icobilist"]
            )

            coxxcar_path = file_paths.get("cobicar")
            icoxxlist_path = file_paths.get("icobilist")
            are_cobis = True
            are_coops = False

        elif self.fingerprint_for.upper() == "COOP":
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "coopcar", "icooplist"]
            )

            coxxcar_path = file_paths.get("coopcar")
            icoxxlist_path = file_paths.get("icooplist")
            are_cobis = False
            are_coops = True

        else:
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "cohpcar", "icohplist"]
            )

            coxxcar_path = file_paths.get("cohpcar")
            icoxxlist_path = file_paths.get("icohplist")
            are_cobis = False
            are_coops = False

        coxx = FeaturizeCOXX(
            path_to_coxxcar=str(coxxcar_path),
            path_to_icoxxlist=str(icoxxlist_path),
            path_to_structure=str(file_paths.get("structure")),
            feature_type=self.feature_type,
            e_range=self.e_range,
            are_coops=are_coops,
            are_cobis=are_cobis,
        )

        return coxx.get_coxx_fingerprint_df(
            spin_type=self.spin_type,
            n_bins=self.n_bins,
            normalize=self.normalize,
            label_list=self.label_list,
        )

    def _get_fingerprints_df(self) -> pd.DataFrame:
        """
        Generate fingerprint objects dataframe using BatchCoxxFingerprint._fingerprint_df.

        Returns:
            A pandas dataframe with COXX fingerprint objects

        """
        paths = [
            os.path.join(self.path_to_lobster_calcs, f)
            for f in os.listdir(self.path_to_lobster_calcs)
            if not f.startswith("t")
            and not f.startswith(".")
            and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
        ]

        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(
                total=len(paths),
                desc=f"Generating {self.fingerprint_for.upper()} fingerprints",
            ) as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._fingerprint_df, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df = pd.concat(row)
        df.sort_index(inplace=True)  # noqa: PD002

        return df


class BatchStructureGraphs:
    """
    Batch Featurizer that generates structure graphs with lobster data.

    :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
    :param add_additional_data_sg: bool indicating whether to include `icoop` and `icobi` data as edge properties
    :param which_bonds: selects which kind of bonds are analyzed. "all" is the default
    :param cutoff_icohp: only bonds that are stronger than cutoff_icohp * strongest ICOHP will be considered.
    :param noise_cutoff: if provided hardcodes the lower limit of icohps considered.
    :param start: start energy for bonding antibonding percent integration
    :param n_jobs: parallel processes to run

    """

    def __init__(
        self,
        path_to_lobster_calcs: str | Path,
        add_additional_data_sg: bool = True,
        which_bonds: str = "all",
        cutoff_icohp: float = 0.10,
        noise_cutoff: float = 0.1,
        start: float | None = None,
        n_jobs: int = 4,
    ):
        """
        Generate structure graphs with LOBSTER data via multiprocessing.

        :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
        :param add_additional_data_sg: bool indicating whether to include `icoop` and `icobi` data as edge properties
        :param which_bonds: selects which kind of bonds are analyzed. "all" is the default
        :param cutoff_icohp: only bonds that are stronger than cutoff_icohp * strongest ICOHP will be considered.
        :param noise_cutoff: if provided hardcodes the lower limit of icohps considered.
        :param start: start energy for bonding antibonding percent integration
        :param n_jobs: parallel processes to run

        """
        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.add_additional_data_sg = add_additional_data_sg
        self.which_bonds = which_bonds
        self.cutff_icohp = cutoff_icohp
        self.noise_cutoff = noise_cutoff
        self.start = start
        self.n_jobs = n_jobs

    def _get_sg_df(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Generate a structure graph with LOBSTER data bonding analysis data.

        :param path_to_lobster_calc: path to root LOBSTER calculation directory

        Returns:
            A  structure graph with LOBSTER data as edge and node properties in structure graph objects
        """
        dir_name = Path(path_to_lobster_calc)
        file_paths = get_file_paths(
            path_to_lobster_calc=path_to_lobster_calc,
            requested_files=["charge", "cohpcar", "icohplist", "icooplist", "icobilist", "madelung", "structure"],
        )

        graph = LobsterGraph(
            path_to_poscar=str(file_paths.get("structure")),
            path_to_charge=str(file_paths.get("charge")),
            path_to_cohpcar=str(file_paths.get("cohpcar")),
            path_to_icohplist=str(file_paths.get("icohplist")),
            add_additional_data_sg=self.add_additional_data_sg,
            path_to_icooplist=str(file_paths.get("icooplist")),
            path_to_icobilist=str(file_paths.get("icobilist")),
            path_to_madelung=str(file_paths.get("madelung")),
            which_bonds=self.which_bonds,
            cutoff_icohp=self.cutff_icohp,
            noise_cutoff=self.noise_cutoff,
            start=self.start,
        )

        ids = dir_name.name

        df = pd.DataFrame(index=[ids])

        df.loc[ids, "structure_graph"] = graph.sg

        return df

    def get_df(self) -> pd.DataFrame:
        """
        Generate a pandas dataframe with structure graph with LOBSTER data.

        Uses multiprocessing to speed up the process.

        Returns:
            Returns a pandas dataframe

        """
        paths = [
            os.path.join(self.path_to_lobster_calcs, f)
            for f in os.listdir(self.path_to_lobster_calcs)
            if not f.startswith("t")
            and not f.startswith(".")
            and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
        ]
        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(paths), desc="Generating Structure Graphs") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._get_sg_df, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df_sg = pd.concat(row)
        df_sg.sort_index(inplace=True)  # noqa: PD002

        return df_sg


class BatchDosFeaturizer:
    """
    BatchFeaturizer to generate Lobster DOS moment features and fingerprints.

    :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
    :param add_element_dos_moments: add element dos moment features alongside orbital dos
    :param normalize: bool to state to normalize the fingerprint data
    :param n_bins: sets number for bins for fingerprint objects
    :param e_range: range of energy relative to fermi for which moment features needs to be computed
    :param n_jobs: number of parallel processes to run
    :param fingerprint_type: Specify fingerprint type to compute, can accept `{s/p/d/f/}summed_{pdos/tdos}`
        (default is summed_pdos)
    :param use_lso_dos: Will force featurizer to use DOSCAR.LSO.lobster instead of DOSCAR.lobster

    """

    def __init__(
        self,
        path_to_lobster_calcs: str | Path,
        add_element_dos_moments: bool = False,
        fingerprint_type: str = "summed_pdos",
        normalize: bool = True,
        n_bins: int = 56,
        e_range: list[float] = [-15.0, 0.0],
        n_jobs=4,
        use_lso_dos: bool = True,
    ):
        """
        Initialize BatchDosFeaturizer.

        :param path_to_lobster_calcs: path to root directory consisting of all lobster calc
        :param add_element_dos_moments: add element dos moment features alongside orbital dos
        :param normalize: bool to state to normalize the fingerprint data
        :param n_bins: sets number for bins for fingerprint objects
        :param e_range: range of energy relative to fermi for which moment features needs to be computed
        :param n_jobs: number of parallel processes to run
        :param fingerprint_type: Specify fingerprint type to compute, can accept `{s/p/d/f/}summed_{pdos/tdos}`
            (default is summed_pdos)
        :param use_lso_dos: Will force featurizer to use DOSCAR.LSO.lobster instead of DOSCAR.lobster
        """
        self.path_to_lobster_calcs = path_to_lobster_calcs
        self.add_element_dos_moments = add_element_dos_moments
        self.fingerprint_type = fingerprint_type
        self.e_range = e_range
        self.normalize = normalize
        self.n_jobs = n_jobs
        self.n_bins = n_bins
        self.use_lso_dos = use_lso_dos

    def _get_dos_moments_df(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Featurize DOSCAR.lobster data using FeaturizeDOSCAR.

        Returns:
            A pandas dataframe with computed PDOS moment features
        """
        file_paths = get_file_paths(
            path_to_lobster_calc=path_to_lobster_calc,
            requested_files=["structure", "doscar"],
            use_lso_dos=self.use_lso_dos,
        )

        featurize_dos = FeaturizeDoscar(
            path_to_doscar=str(file_paths.get("doscar")),
            path_to_structure=str(file_paths.get("structure")),
            add_element_dos_moments=self.add_element_dos_moments,
            e_range=self.e_range,
        )

        return featurize_dos.get_df()

    def _get_dos_fingerprints_df(self, path_to_lobster_calc: str | Path) -> pd.DataFrame:
        """
        Featurize DOSCAR.lobster data into fingerprints using FeaturizeDOSCAR.

        :param path_to_lobster_calc: path to root LOBSTER calculation directory.

        Returns:
            A pandas dataframe with DOS fingerprint objects
        """
        file_paths = get_file_paths(
            path_to_lobster_calc=path_to_lobster_calc,
            requested_files=["structure", "doscar"],
            use_lso_dos=self.use_lso_dos,
        )

        featurize_dos = FeaturizeDoscar(
            path_to_doscar=str(file_paths.get("doscar")),
            path_to_structure=str(file_paths.get("structure")),
            e_range=self.e_range,
        )

        return featurize_dos.get_fingerprint_df(
            fp_type=self.fingerprint_type,
            normalize=self.normalize,
            n_bins=self.n_bins,
        )

    def get_df(self) -> pd.DataFrame:
        """
        Generate a pandas dataframe with all moment features.

        Moment features are PDOS (optional: element dos) center, width, skewness, kurtosis
        and upper band edge.

        Returns:
            A pandas dataframe with moment features
        """
        paths = [
            os.path.join(self.path_to_lobster_calcs, f)
            for f in os.listdir(self.path_to_lobster_calcs)
            if not f.startswith("t")
            and not f.startswith(".")
            and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
        ]
        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(paths), desc="Generating PDOS moment features") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._get_dos_moments_df, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df_dos = pd.concat(row)
        df_dos.sort_index(inplace=True)  # noqa: PD002

        return df_dos

    def get_fingerprints_df(self) -> pd.DataFrame:
        """
        Generate a pandas dataframe with DOS fingerprints.

        Returns:
            A pandas dataframe with fingerprint objects
        """
        paths = [
            os.path.join(self.path_to_lobster_calcs, f)
            for f in os.listdir(self.path_to_lobster_calcs)
            if not f.startswith("t")
            and not f.startswith(".")
            and os.path.isdir(os.path.join(self.path_to_lobster_calcs, f))
        ]
        row = []
        with (
            mp.Pool(processes=self.n_jobs, maxtasksperchild=1) as pool,
            tqdm(total=len(paths), desc="Generating DOS fingerprints") as pbar,
        ):
            for _, result in enumerate(pool.imap_unordered(self._get_dos_fingerprints_df, paths, chunksize=1)):
                pbar.update()
                row.append(result)

        df_dos_fp = pd.concat(row)
        df_dos_fp.sort_index(inplace=True)  # noqa: PD002

        return df_dos_fp
