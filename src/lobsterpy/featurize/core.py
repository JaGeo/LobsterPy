# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module defines classes to featurize Lobster data ready to be used for ML studies."""

from __future__ import annotations

import gzip
import json
import warnings
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

try:
    from mendeleev import element
except ImportError:
    element = None
from monty.dev import requires
from numpy import ndarray
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster import Charge, Doscar, Grosspop, Icohplist, MadelungEnergies
from scipy.integrate import trapezoid
from scipy.signal import hilbert
from scipy.stats import kurtosis, skew, wasserstein_distance

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.featurize.utils import CoxxFingerprint, get_file_paths

warnings.filterwarnings("ignore")


class FeaturizeLobsterpy:
    """
    Class to featurize lobsterpy data.

    :param path_to_lobster_calc: path to parent directory containing lobster calc outputs
    :param path_to_json: path to lobster lightweight json
    :param bonds: "all" or "cation-anion" bonds
    :param orbital_resolved: bool indicating whether LobsterPy analysis is performed orbital wise
    :param analysis_kwargs: optional keyword arguments to be passed to the Analysis class
    """

    def __init__(
        self,
        path_to_lobster_calc: str | Path | None = None,
        path_to_json: str | Path | None = None,
        orbital_resolved: bool = False,
        bonds: Literal["cation-anion", "all"] = "all",
        **analysis_kwargs,
    ):
        """Initialize featurizer."""
        self.path_to_json = path_to_json
        self.path_to_lobster_calc = path_to_lobster_calc
        self.orbital_resolved = orbital_resolved
        self.bonds = bonds
        self.analysis_kwargs = analysis_kwargs

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        Featurize LobsterPy condensed bonding analysis data.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.

        Returns:
            Returns a pandas dataframe with lobsterpy icohp statistics
        """
        if self.path_to_json and not self.path_to_lobster_calc:
            # read the lightweight lobster json files using read_lobster_lightweight_json method
            data = FeaturizeLobsterpy.read_lobster_lightweight_json(path_to_json=self.path_to_json)
            type_pop = "COHP"
            if not ids:
                ids = Path(self.path_to_json).name.split(".")[0]

        elif self.path_to_lobster_calc and not self.path_to_json:
            # get lobsterpy condensed bonding analysis data using get_lobsterpy_cba_dict method
            type_pop, data = FeaturizeLobsterpy.get_lobsterpy_cba_dict(
                path_to_lobster_calc=self.path_to_lobster_calc,
                bonds=self.bonds,
                orbital_resolved=self.orbital_resolved,
                **self.analysis_kwargs,
            )

            if not ids:
                ids = Path(self.path_to_lobster_calc).name

        else:
            raise ValueError("Please provide either path to lightweight lobster jsons or path to lobster calc")
        # define a pandas dataframe
        df = pd.DataFrame(index=[ids])

        icohp_mean = []
        icohp_sum = []
        icohp_mean_orb_bndg = []
        icohp_sum_orb_bndg = []
        icohp_mean_orb_antibndg = []
        icohp_sum_orb_antibndg = []
        bond_orb = []
        antibond_orb = []
        bond = []
        antibond = []

        # extract lobsterpy icohp related data for bond type specified

        # Results will differ for "all" and "cation-anion" mode.

        # In "all" bonds mode, the bonds will come up twice, also
        # cation-cation, anion-anion bonds will also be considered

        if self.bonds == "all":
            bond_type = "all_bonds"
        elif self.bonds == "cation-anion":
            bond_type = "cation_anion_bonds"

        if bond_type == "cation_anion_bonds" and (
            not data[bond_type]
            or not data[bond_type]["lobsterpy_data"]
            or not data[bond_type]["lobsterpy_data"]["sites"]
        ):
            raise Exception(f"No {self.bonds} bonds detected for {ids} structure. Please switch to `all` bonds mode")

        if data[bond_type]["lobsterpy_data"]:
            for site_data in data[bond_type]["lobsterpy_data"]["sites"].values():
                if site_data["bonds"]:
                    for bond_data in site_data["bonds"].values():
                        icohp_mean.append(float(bond_data[f"I{type_pop}_mean"]))
                        icohp_sum.append(float(bond_data[f"I{type_pop}_sum"]))
                        bond.append(bond_data["bonding"]["perc"])
                        antibond.append(bond_data["antibonding"]["perc"])
                        if self.orbital_resolved:
                            if (
                                bond_data["orbital_data"]["orbital_summary_stats"]
                                and "max_bonding_contribution" in bond_data["orbital_data"]["orbital_summary_stats"]
                            ):
                                for orb_pair in bond_data["orbital_data"]["orbital_summary_stats"][
                                    "max_bonding_contribution"
                                ]:
                                    icohp_mean_orb_bndg.append(bond_data["orbital_data"][orb_pair][f"I{type_pop}_mean"])
                                    icohp_sum_orb_bndg.append(bond_data["orbital_data"][orb_pair][f"I{type_pop}_sum"])
                                    bond_orb.append(
                                        bond_data["orbital_data"][orb_pair]["orb_contribution_perc_bonding"]
                                    )
                            if (
                                bond_data["orbital_data"]["orbital_summary_stats"]
                                and "max_antibonding_contribution" in bond_data["orbital_data"]["orbital_summary_stats"]
                            ):
                                for orb_pair in bond_data["orbital_data"]["orbital_summary_stats"][
                                    "max_antibonding_contribution"
                                ]:
                                    icohp_mean_orb_antibndg.append(
                                        bond_data["orbital_data"][orb_pair][f"I{type_pop}_mean"]
                                    )
                                    icohp_sum_orb_antibndg.append(
                                        bond_data["orbital_data"][orb_pair][f"I{type_pop}_sum"]
                                    )
                                    antibond_orb.append(
                                        bond_data["orbital_data"][orb_pair]["orb_contribution_perc_antibonding"]
                                    )

        # add ICOHP stats data (mean, min, max, standard deviation) as columns to the dataframe

        df.loc[ids, f"I{type_pop}_mean_avg"] = 0 if len(icohp_mean) == 0 else np.mean(icohp_mean)
        df.loc[ids, f"I{type_pop}_mean_max"] = 0 if len(icohp_mean) == 0 else np.max(icohp_mean)
        df.loc[ids, f"I{type_pop}_mean_min"] = 0 if len(icohp_mean) == 0 else np.min(icohp_mean)
        df.loc[ids, f"I{type_pop}_mean_std"] = 0 if len(icohp_mean) == 0 else np.std(icohp_mean)

        df.loc[ids, f"I{type_pop}_sum_avg"] = 0 if len(icohp_sum) == 0 else np.mean(icohp_sum)
        df.loc[ids, f"I{type_pop}_sum_max"] = 0 if len(icohp_sum) == 0 else np.max(icohp_sum)
        df.loc[ids, f"I{type_pop}_sum_min"] = 0 if len(icohp_sum) == 0 else np.min(icohp_sum)
        df.loc[ids, f"I{type_pop}_sum_std"] = 0 if len(icohp_sum) == 0 else np.std(icohp_sum)

        df.loc[ids, "bonding_perc_avg"] = 0 if len(bond) == 0 else np.mean(bond)
        df.loc[ids, "bonding_perc_max"] = 0 if len(bond) == 0 else np.max(bond)
        df.loc[ids, "bonding_perc_min"] = 0 if len(bond) == 0 else np.min(bond)
        df.loc[ids, "bonding_perc_std"] = 0 if len(bond) == 0 else np.std(bond)

        df.loc[ids, "antibonding_perc_avg"] = 0 if len(antibond) == 0 else np.mean(antibond)
        df.loc[ids, "antibonding_perc_min"] = 0 if len(antibond) == 0 else np.min(antibond)
        df.loc[ids, "antibonding_perc_max"] = 0 if len(antibond) == 0 else np.max(antibond)
        df.loc[ids, "antibonding_perc_std"] = 0 if len(antibond) == 0 else np.std(antibond)

        if self.orbital_resolved:
            # bonding orbital
            df.loc[ids, f"I{type_pop}_bndg_orb_mean_avg"] = (
                0 if len(icohp_mean_orb_bndg) == 0 else np.mean(icohp_mean_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_mean_max"] = (
                0 if len(icohp_mean_orb_bndg) == 0 else np.max(icohp_mean_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_mean_min"] = (
                0 if len(icohp_mean_orb_bndg) == 0 else np.min(icohp_mean_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_mean_std"] = (
                0 if len(icohp_mean_orb_bndg) == 0 else np.std(icohp_mean_orb_bndg)
            )

            df.loc[ids, f"I{type_pop}_bndg_orb_sum_avg"] = (
                0 if len(icohp_sum_orb_bndg) == 0 else np.mean(icohp_sum_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_sum_max"] = (
                0 if len(icohp_sum_orb_bndg) == 0 else np.max(icohp_sum_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_sum_min"] = (
                0 if len(icohp_sum_orb_bndg) == 0 else np.min(icohp_sum_orb_bndg)
            )
            df.loc[ids, f"I{type_pop}_bndg_orb_sum_std"] = (
                0 if len(icohp_sum_orb_bndg) == 0 else np.std(icohp_sum_orb_bndg)
            )

            df.loc[ids, "bonding_orb_perc_avg"] = 0 if len(bond_orb) == 0 else np.mean(bond_orb)
            df.loc[ids, "bonding_orb_perc_max"] = 0 if len(bond_orb) == 0 else np.max(bond_orb)
            df.loc[ids, "bonding_orb_perc_min"] = 0 if len(bond_orb) == 0 else np.min(bond_orb)
            df.loc[ids, "bonding_orb_perc_std"] = 0 if len(bond_orb) == 0 else np.std(bond_orb)

            # anti-bonding orbital
            df.loc[ids, f"I{type_pop}_antibndg_orb_mean_avg"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.mean(icohp_mean_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_mean_max"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.max(icohp_mean_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_mean_min"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.min(icohp_mean_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_mean_std"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.std(icohp_mean_orb_antibndg)
            )

            df.loc[ids, f"I{type_pop}_antibndg_orb_sum_avg"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.mean(icohp_sum_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_sum_max"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.max(icohp_sum_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_sum_min"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.min(icohp_sum_orb_antibndg)
            )
            df.loc[ids, f"I{type_pop}_antibndg_orb_sum_std"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.std(icohp_sum_orb_antibndg)
            )

            df.loc[ids, "antibonding_orb_perc_avg"] = 0 if len(antibond_orb) == 0 else np.mean(antibond_orb)
            df.loc[ids, "antibonding_orb_perc_max"] = 0 if len(antibond_orb) == 0 else np.max(antibond_orb)
            df.loc[ids, "antibonding_orb_perc_min"] = 0 if len(antibond_orb) == 0 else np.min(antibond_orb)
            df.loc[ids, "antibonding_orb_perc_std"] = 0 if len(antibond_orb) == 0 else np.std(antibond_orb)

        # add madelung energies for the structure
        df.loc[ids, "Madelung_Mull"] = data["madelung_energies"]["Mulliken"]
        df.loc[ids, "Madelung_Loew"] = data["madelung_energies"]["Loewdin"]

        return df

    @staticmethod
    def read_lobster_lightweight_json(path_to_json: str | Path) -> dict:
        """
        Read the lightweight JSON.gz files and return a Python dictionary object.

        :param path_to_json: path to lobsterpy lightweight json file

        Returns:
            Returns a dictionary with lobster summarized bonding analysis data

        """
        with gzip.open(str(path_to_json), "rb") as f:
            data = json.loads(f.read().decode("utf-8"))

        lobster_data = {}
        for item in data:
            for key in item:
                # check applicable only for updated cba jsons from atomate2
                try:
                    if (key == "cation_anion_bonds" or key == "all_bonds") and (
                        "sites" in item[key]["lobsterpy_data"]["sites"]
                    ):
                        item[key]["lobsterpy_data"]["sites"] = item[key]["lobsterpy_data"]["sites"]["sites"]
                except KeyError:
                    pass
            lobster_data.update(item)

        return lobster_data

    @staticmethod
    def get_lobsterpy_cba_dict(
        path_to_lobster_calc: str | Path, bonds: str, orbital_resolved: bool, **analysis_kwargs
    ) -> tuple[str, dict]:
        """
        Generate a Python dictionary object using the Analysis class with condensed bonding analysis data.

        :param path_to_lobster_calc: path to lobsterpy lightweight json file
        :param bonds: "all" or "cation-anion" bonds
        :param orbital_resolved:  bool indicating whether analysis is performed orbital wise
        :param analysis_kwargs: optional keyword arguments to be passed to the Analysis class

        Returns:
            Returns a string indicating type of population analyzed and a
            dictionary with lobster summarized bonding analysis data

        """
        which_bonds = bonds.replace("-", "_")
        bond_type = f"{which_bonds}_bonds"

        # Initialize a set of default kwargs of Analysis class
        default_kwargs = {
            "are_cobis": analysis_kwargs.get("are_cobis", False),
            "are_coops": analysis_kwargs.get("are_coops", False),
            "type_charge": analysis_kwargs.get("type_charge", "Mulliken"),
            "cutoff_icohp": analysis_kwargs.get("cutoff_icohp", 0.10),
            "noise_cutoff": analysis_kwargs.get("noise_cutoff", 0.1),
            "orbital_resolved": orbital_resolved,
            "start": analysis_kwargs.get("start"),
            "summed_spins": analysis_kwargs.get("summed_spins", False),
            "which_bonds": which_bonds,
        }

        if default_kwargs.get("are_cobis", False):
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc,
                requested_files=["structure", "cobicar", "icobilist", "charge"],
            )
            default_kwargs.update(
                {
                    "path_to_poscar": str(file_paths.get("structure")),
                    "path_to_icohplist": str(file_paths.get("icobilist")),
                    "path_to_cohpcar": str(file_paths.get("cobicar")),
                    "path_to_charge": str(file_paths.get("charge")),
                }
            )

        elif default_kwargs.get("are_coops", False):
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc,
                requested_files=["structure", "coopcar", "icooplist", "charge"],
            )
            default_kwargs.update(
                {
                    "path_to_poscar": str(file_paths.get("structure")),
                    "path_to_icohplist": str(file_paths.get("icooplist")),
                    "path_to_cohpcar": str(file_paths.get("coopcar")),
                    "path_to_charge": str(file_paths.get("charge")),
                }
            )
        else:
            file_paths = get_file_paths(
                path_to_lobster_calc=path_to_lobster_calc,
                requested_files=["structure", "cohpcar", "icohplist", "charge"],
            )
            default_kwargs.update(
                {
                    "path_to_poscar": str(file_paths.get("structure")),
                    "path_to_icohplist": str(file_paths.get("icohplist")),
                    "path_to_cohpcar": str(file_paths.get("cohpcar")),
                    "path_to_charge": str(file_paths.get("charge")),
                }
            )

        try:
            analyse = Analysis(**default_kwargs)
            type_pop = analyse._get_pop_type()

            data = {bond_type: {"lobsterpy_data": analyse.condensed_bonding_analysis}}
        except ValueError:
            data = {bond_type: {"lobsterpy_data": {}}}

            if default_kwargs.get("are_cobis"):
                type_pop = "COBI"
            elif default_kwargs.get("are_coops"):
                type_pop = "COOP"
            else:
                type_pop = "COHP"

        try:
            madelung_path = get_file_paths(path_to_lobster_calc=path_to_lobster_calc, requested_files=["madelung"])
            madelung_obj = MadelungEnergies(filename=str(madelung_path.get("madelung")))

            madelung_energies = {
                "Mulliken": madelung_obj.madelungenergies_Mulliken,
                "Loewdin": madelung_obj.madelungenergies_Loewdin,
                "Ewald_splitting": madelung_obj.ewald_splitting,
            }
            data["madelung_energies"] = madelung_energies
        except Exception:
            warnings.warn(
                "MadelungEnergies.lobster file not found in Lobster calc directory provided"
                " Will set Madelung Energies for crystal structure values to NaN",
                stacklevel=2,
            )
            madelung_energies = {
                "Mulliken": np.nan,
                "Loewdin": np.nan,
                "Ewald_splitting": np.nan,
            }

            data["madelung_energies"] = madelung_energies

        return type_pop, data


class FeaturizeCOXX:
    """
    Class to featurize COHPCAR, COBICAR or COOPCAR data.

    :param path_to_coxxcar: path to COXXCAR.lobster (e.g., `COXXCAR.lobster`)
    :param path_to_icoxxlist: path to ICOXXLIST.lobster (e.g., `ICOXXLIST.lobster`)
    :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
    :param feature_type: set the feature type for moment features and fingerprints.
        Possible options are `bonding`, `antibonding` or `overall`.
    :param are_cobis: bool indicating if file contains COBI/ICOBI data.
    :param are_coops: bool indicating if file contains COOP/ICOOP data.
    :param e_range: range of energy relative to fermi for which moment features needs to be computed
    """

    def __init__(
        self,
        path_to_coxxcar: str | Path,
        path_to_icoxxlist: str | Path,
        path_to_structure: str | Path,
        feature_type: Literal["bonding", "antibonding", "overall"],
        e_range: list[float] = [-10.0, 0.0],
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
        """
        Featurize COHPCAR, COBICAR or COOPCAR data.

        :param path_to_coxxcar: path to COXXCAR.lobster (e.g., `COXXCAR.lobster`)
        :param path_to_icoxxlist: path to ICOXXLIST.lobster (e.g., `ICOXXLIST.lobster`)
        :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
        :param feature_type: set the feature type for moment features and fingerprints.
            Possible options are `bonding`, `antibonding` or `overall`.
        :param are_cobis: bool indicating if file contains COBI/ICOBI data
        :param are_coops: bool indicating if file contains COOP/ICOOP data
        :param e_range: range of energy relative to fermi for which moment features needs to be computed

        """
        self.path_to_coxxcar = path_to_coxxcar
        self.path_to_icoxxlist = path_to_icoxxlist
        self.path_to_structure = path_to_structure
        self.feature_type = feature_type
        self.e_range = e_range
        self.are_cobis = are_cobis
        self.are_coops = are_coops
        self.icoxxlist = Icohplist(
            filename=self.path_to_icoxxlist,
            are_cobis=self.are_cobis,
            are_coops=self.are_coops,
        )
        self.completecoxx = CompleteCohp.from_file(
            fmt="LOBSTER",
            filename=self.path_to_coxxcar,
            structure_file=self.path_to_structure,
            are_cobis=self.are_cobis,
            are_coops=self.are_coops,
        )

    def get_coxx_fingerprint_df(
        self,
        ids: str | None = None,
        label_list: list[str] | None = None,
        orbital: str | None = None,
        per_bond: bool = True,
        spin_type: Literal["up", "down", "summed"] = "summed",
        binning: bool = True,
        n_bins: int = 56,
        normalize: bool = True,
    ) -> pd.DataFrame:
        """
        Generate a COXX fingerprints dataframe.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.
        :param spin_type: Specify spin type. Can accept '{summed/up/down}'
            (default is summed)
        :param binning: If true coxxs will be binned
        :param n_bins: Number of bins to be used in the fingerprint (default is 256)
        :param normalize: If true, normalizes the area under fp to equal to 1 (default is True)
        :param label_list: Specify bond labels as a list for which cohp fingerprints are needed
        :param orbital: Orbital for which fingerprint needs is to be computed. Cannot be used independently.
            Always a needs label_list.
        :param per_bond: Will scale cohp values by number of bonds i.e. length of label_list arg
            (Only affects when label_list is not None)

        Raises:
            Exception: If spin_type is not one of the accepted values {summed/up/down}.
            ValueError: If feature_type is not one of the accepted values {bonding/antibonding/overall}.

        Returns:
            A pandas dataframe with the COXX fingerprint
            of format (energies, coxx, fp_type, spin_type, n_bins, bin_width)

        """
        coxxcar_obj = self.completecoxx

        energies = coxxcar_obj.energies - coxxcar_obj.efermi

        min_e = self.e_range[0]
        max_e = self.e_range[-1]

        if max_e is None:
            max_e = np.max(energies)

        if min_e is None:
            min_e = np.min(energies)

        if label_list is not None:
            divisor = len(label_list) if per_bond else 1
        if label_list is not None and orbital is None:
            coxxcar_obj = coxxcar_obj.get_summed_cohp_by_label_list(label_list, divisor=divisor).get_cohp()
        elif label_list and orbital:
            coxxcar_obj = coxxcar_obj.get_summed_cohp_by_label_and_orbital_list(
                label_list=label_list,
                orbital_list=[orbital] * len(label_list),
                divisor=divisor,
            ).get_cohp()
        else:
            coxxcar_obj = coxxcar_obj.get_cohp()

        if spin_type == "up":
            coxx_all = coxxcar_obj[Spin.up]
        elif spin_type == "down":
            if Spin.down in coxxcar_obj:
                coxx_all = coxxcar_obj[Spin.down]
            else:
                raise ValueError("LOBSTER calculation is non-spin polarized. Please switch spin_type to `up`")
        elif spin_type == "summed":
            if Spin.down in coxxcar_obj:
                coxx_all = coxxcar_obj[Spin.up] + coxxcar_obj[Spin.down]
            else:
                coxx_all = coxxcar_obj[Spin.up]
        else:
            raise Exception("Check the spin_type argument.Possible options are summed/up/down")

        coxx_dict = {}
        tol = 1e-6
        if not self.are_cobis and not self.are_coops:
            coxx_dict["bonding"] = np.array([scohp if scohp <= tol else 0 for scohp in coxx_all])
            coxx_dict["antibonding"] = np.array([scohp if scohp >= tol else 0 for scohp in coxx_all])
        else:
            coxx_dict["antibonding"] = np.array([scohp if scohp <= tol else 0 for scohp in coxx_all])
            coxx_dict["bonding"] = np.array([scohp if scohp >= tol else 0 for scohp in coxx_all])

        coxx_dict["overall"] = coxx_all

        try:
            if ids:
                df = pd.DataFrame(index=[ids], columns=["COXX_FP"])
            else:
                ids = Path(self.path_to_coxxcar).parent.name
                df = pd.DataFrame(index=[ids], columns=["COXX_FP"])

            coxxs = coxx_dict[self.feature_type]
            if len(energies) < n_bins:
                inds = np.where((energies >= min_e - tol) & (energies <= max_e + tol))
                fp = CoxxFingerprint(
                    energies[inds],
                    coxxs[inds],
                    self.feature_type,
                    spin_type,
                    len(energies),
                    np.diff(energies)[0],
                )

                df.loc[ids, "COXX_FP"] = fp

                return df

            if binning:
                ener_bounds = np.linspace(min_e, max_e, n_bins + 1)
                ener = ener_bounds[:-1] + (ener_bounds[1] - ener_bounds[0]) / 2.0
                bin_width = np.diff(ener)[0]
            else:
                ener_bounds = np.array(energies)
                ener = np.append(energies, [energies[-1] + np.abs(energies[-1]) / 10])
                n_bins = len(energies)
                bin_width = np.diff(energies)[0]

            coxx_rebin = np.zeros(ener.shape)

            for ii, e1, e2 in zip(range(len(ener)), ener_bounds[0:-1], ener_bounds[1:]):
                inds = np.where((energies >= e1) & (energies < e2))
                coxx_rebin[ii] = np.sum(coxxs[inds])
            if normalize:  # scale DOS bins to make area under histogram equal 1
                area = np.sum([abs(coxx) * bin_width for coxx in coxx_rebin])
                coxx_rebin_sc = coxx_rebin / area
            else:
                coxx_rebin_sc = coxx_rebin

            fp = CoxxFingerprint(
                np.array([ener]),
                coxx_rebin_sc,
                self.feature_type,
                spin_type,
                n_bins,
                bin_width,
            )

            df.loc[ids, "COXX_FP"] = fp

            return df

        except KeyError:
            raise ValueError(
                "Please recheck fingerprint type requested argument. Possible options are bonding/antibonding/overall"
            )

    def _calculate_wicoxx_ein(self) -> tuple[float, float, float, float]:
        r"""
        Calculate weighted ICOXX (xx ==hp,op,bi) and EIN of crystal structure based on work by Nelson et al.

        More details can be found here: https://doi.org/10.1016/B978-0-12-823144-9.00120-5

        \\[\bar{COHP}(E) = \\sum_{i} \\left( w_{i} \\cdot COHP_{i}(E) \right)\\]

        where:
        - \\(w_{i} = \frac{ICOHP_{i}}{ICOHP_{\text{total}}}\\)

        \\[EIN = \frac{{ICOHP_{\text{{total}}}}}{{\\overline{ICOHP}}} \\cdot \frac{2}{{N_{\text{{Atoms}}}}}\\]

        where:
        - \\(\\overline{ICOHP}\\) represents the average ICOHP value.
        - \\(ICOHP_{\text{{total}}}\\) is the total ICOHP value.
        - \\(N_{\text{{Atoms}}}\\) is the number of atoms in the system.

        Returns:
            Percent bonding, Percent anti-bonding, weighted icoxx, effective interaction number

        """
        list_labels = list(self.icoxxlist.icohplist.keys())
        # Compute sum of icohps
        icoxx_total = self.icoxxlist.icohpcollection.get_summed_icohp_by_label_list(list_labels)

        summed_weighted_coxx = []

        for lab in list_labels:
            for cohp in self.completecoxx.get_cohp_by_label(f"{lab}", summed_spin_channels=True).cohp.values():
                coxx = cohp
            icoxx = self.icoxxlist.icohpcollection.get_icohp_by_label(lab, summed_spin_channels=True)
            # calculate the weights based on icohp contri to total icohp of the structure
            try:
                weight = icoxx / icoxx_total
            except ZeroDivisionError:
                weight = 0
            weighted_coxx = weight * coxx
            summed_weighted_coxx.append(weighted_coxx)

        summed = np.sum(summed_weighted_coxx, axis=0)
        # below fermi
        indices = self.completecoxx.energies <= self.completecoxx.efermi
        en_bf = self.completecoxx.energies[indices]
        coxx_bf = summed[indices]

        w_icoxx = trapezoid(coxx_bf, en_bf)

        ein = (icoxx_total / w_icoxx) * (2 / self.completecoxx.structure.num_sites)  # calc effective interaction number

        # percent bonding-anitbonding
        tol = 1e-6
        if not self.icoxxlist.are_cobis and not self.icoxxlist.are_coops:
            bonding_indices = coxx_bf <= tol
            antibonding_indices = coxx_bf >= tol
            bnd = abs(trapezoid(en_bf[bonding_indices], coxx_bf[bonding_indices]))
            antibnd = abs(trapezoid(en_bf[antibonding_indices], coxx_bf[antibonding_indices]))
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100

        elif self.icoxxlist.are_cobis or self.icoxxlist.are_coops:
            bonding_indices = coxx_bf >= tol
            antibonding_indices = coxx_bf <= tol
            bnd = abs(trapezoid(coxx_bf[bonding_indices], en_bf[bonding_indices]))
            antibnd = abs(trapezoid(coxx_bf[antibonding_indices], en_bf[antibonding_indices]))
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100

        return per_bnd, per_antibnd, w_icoxx, ein

    @staticmethod
    def _calc_moment_features(
        complete_coxx_obj: CompleteCohp,
        feature_type: Literal["bonding", "antibonding", "overall"],
        e_range: list[float],
        label_list: list[str] | None = None,
        orbital: str | None = None,
        per_bond: bool = True,
    ) -> tuple[float, float, float, float, float]:
        """
        Calculate band center,width, skewness, and kurtosis of the COXX.

        :param complete_coxx_obj: CompleteCohp object
        :param feature_type: feature type for moment features calculation
        :param e_range: range of energy relative to fermi for which moment features needs to be computed
        :param label_list: List of bond labels
        :param orbital: orbital for which moment features need to be calculated. Cannot be used independently.
            Always needs a label_list.
        :param per_bond: Will scale cohp values by number of bonds i.e. length of label_list arg
            (Only affects when label_list is not None)

        Returns:
            coxx center, width, skewness, kurtosis, and edge in eV
        """
        if label_list is not None:
            divisor = len(label_list) if per_bond else 1

        if label_list and orbital is None:
            coxxcar = complete_coxx_obj.get_summed_cohp_by_label_list(label_list, divisor=divisor).get_cohp()
        elif label_list and orbital:
            coxxcar = complete_coxx_obj.get_summed_cohp_by_label_and_orbital_list(
                label_list=label_list,
                orbital_list=[orbital] * len(label_list),
                divisor=divisor,
                summed_spin_channels=False,
            ).get_cohp()
        else:
            coxxcar = complete_coxx_obj.get_cohp()

        coxx_all = coxxcar[Spin.up] + coxxcar[Spin.down] if Spin.down in coxxcar else coxxcar[Spin.up]
        energies = complete_coxx_obj.energies - complete_coxx_obj.efermi

        coxx_dict = {}
        tol = 1e-6
        if not complete_coxx_obj.are_cobis and not complete_coxx_obj.are_coops:
            coxx_dict["bonding"] = np.array([scohp if scohp <= tol else 0 for scohp in coxx_all])
            coxx_dict["antibonding"] = np.array([scohp if scohp >= tol else 0 for scohp in coxx_all])
            coxx_dict["overall"] = coxx_all
        else:
            coxx_dict["antibonding"] = np.array([scohp if scohp <= tol else 0 for scohp in coxx_all])
            coxx_dict["bonding"] = np.array([scohp if scohp >= tol else 0 for scohp in coxx_all])
            coxx_dict["overall"] = coxx_all

        try:
            coxx_center = FeaturizeCOXX.get_coxx_center(
                coxx=coxx_dict[feature_type],
                energies=energies,
                e_range=e_range,
            )
            coxx_width = np.sqrt(
                FeaturizeCOXX.get_n_moment(
                    n=2,
                    coxx=coxx_dict[feature_type],
                    energies=energies,
                    e_range=e_range,
                )
            )
            coxx_skew = FeaturizeCOXX.get_n_moment(
                n=3,
                coxx=coxx_dict[feature_type],
                energies=energies,
                e_range=e_range,
            ) / FeaturizeCOXX.get_n_moment(
                n=2,
                coxx=coxx_dict[feature_type],
                energies=energies,
                e_range=e_range,
            ) ** (3 / 2)

            coxx_kurt = (
                FeaturizeCOXX.get_n_moment(
                    n=4,
                    coxx=coxx_dict[feature_type],
                    energies=energies,
                    e_range=e_range,
                )
                / FeaturizeCOXX.get_n_moment(
                    n=2,
                    coxx=coxx_dict[feature_type],
                    energies=energies,
                    e_range=e_range,
                )
                ** 2
            )
            coxx_edge = FeaturizeCOXX.get_cohp_edge(coxx=coxx_dict[feature_type], energies=energies, e_range=e_range)
        except KeyError:
            raise ValueError(
                "Please recheck feature type requested argument. Possible options are bonding/antibonding/overall"
            )

        return coxx_center, coxx_width, coxx_skew, coxx_kurt, coxx_edge

    @staticmethod
    def get_coxx_center(
        coxx: ndarray[np.floating],
        energies: ndarray[np.floating],
        e_range: list[float],
    ) -> float:
        """
        Get the bandwidth, defined as the first moment of the COXX.

        :param coxx: COXX array
        :param energies: energies corresponding  COXX
        :param e_range: range of energy to compute coxx center

        Returns:
            coxx center in eV
        """
        return FeaturizeCOXX.get_n_moment(n=1, coxx=coxx, energies=energies, e_range=e_range, center=False)

    @staticmethod
    def get_n_moment(
        n: float,
        coxx: ndarray[np.floating],
        energies: ndarray[np.floating],
        e_range: list[float] | None,
        center: bool = True,
    ) -> float:
        """
        Get the nth moment of COXX.

        :param n: The order for the moment
        :param coxx: COXX array
        :param energies: energies array
        :param e_range: range of energy to compute nth moment
        :param center: Take moments with respect to the COXX center
        :param e_range: range of energy to compute nth moment

        Returns:
            COXX nth moment in eV
        """
        if e_range:
            min_e = e_range[0]
            max_e = e_range[1]
            if min_e is None:
                min_e = min(energies)
            if max_e is None:
                max_e = max(energies)
        else:
            min_e = min(energies)
            max_e = max(energies)

        tol = 1e-6

        mask = (energies >= min_e - tol) & (energies <= max_e + tol)

        coxx = coxx[mask]
        energies = energies[mask]

        if center:
            coxx_center = FeaturizeCOXX.get_coxx_center(coxx=coxx, energies=energies, e_range=[min_e, max_e])
            p = energies - coxx_center
        else:
            p = energies

        return np.trapezoid(p**n * coxx, x=energies) / np.trapezoid(coxx, x=energies)

    @staticmethod
    def get_cohp_edge(
        coxx: ndarray[np.floating],
        energies: ndarray[np.floating],
        e_range: list[float] | None,
    ):
        """
        Get the highest peak position of hilbert transformed COXX.

        :param coxx: COXX array
        :param energies: energies array
        :param e_range: range of energy to coxx edge (max peak position)

        Returns:
            COXX edge (max peak position) in eV
        """
        if e_range:
            min_e = e_range[0]
            max_e = e_range[1]
            if min_e is None:
                min_e = min(energies)
            if max_e is None:
                max_e = max(energies)
        else:
            min_e = min(energies)
            max_e = max(energies)

        tol = 1e-6

        mask = (energies >= min_e - tol) & (energies <= max_e + tol)

        coxx = coxx[mask]
        energies = energies[mask]

        coxx_transformed = np.imag(hilbert(coxx))

        return energies[np.argmax(coxx_transformed)]

    def get_summarized_coxx_df(
        self,
        ids: str | None = None,
        label_list: list[str] | None = None,
        per_bond: bool = True,
    ) -> pd.DataFrame:
        """
        Get a pandas dataframe with COXX features.

        Features consist of weighted ICOXX, effective interaction number and
        moment features of COXX in the selected energy range.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.
        :param label_list: list of bond labels to be used for generating features from
            COHPCAR.lobster or COOPCAR.lobster or COBICAR.lobster. Default is None.
            When None, all bond labels are used.
        :param per_bond: Defaults to True. When True, features are normalized
            by total number of bonds in the structure.

        Returns:
            Returns a pandas dataframe with cohp/cobi/coop related features as per input file

        """
        if ids:
            df = pd.DataFrame(index=[ids])
        else:
            ids = Path(self.path_to_coxxcar).parent.name
            df = pd.DataFrame(index=[ids])

        (
            per_bnd_xx,
            per_antibnd_xx,
            w_icoxx,
            ein_xx,
        ) = self._calculate_wicoxx_ein()

        if self.are_cobis:
            type_pop = "COBI"
        elif self.are_coops:
            type_pop = "COOP"
        else:
            type_pop = "COHP"

        cc, cw, cs, ck, ce = FeaturizeCOXX._calc_moment_features(
            complete_coxx_obj=self.completecoxx,
            label_list=label_list,
            per_bond=per_bond,
            feature_type=self.feature_type,
            e_range=self.e_range,
            orbital=None,
        )
        df.loc[ids, f"bnd_wI{type_pop}"] = per_bnd_xx
        df.loc[ids, f"antibnd_wI{type_pop}"] = per_antibnd_xx
        df.loc[ids, f"w_I{type_pop}"] = w_icoxx
        df.loc[ids, f"EIN_I{type_pop}"] = ein_xx
        df.loc[ids, f"center_{type_pop}"] = cc
        df.loc[ids, f"width_{type_pop}"] = cw
        df.loc[ids, f"skewness_{type_pop}"] = cs
        df.loc[ids, f"kurtosis_{type_pop}"] = ck
        df.loc[ids, f"edge_{type_pop}"] = ce

        return df


@requires(
    element is not None,
    "FeaturizeCharges requires mendeleev. Reinstall package with `pip install lobsterpy[featurizer]`.",
)
class FeaturizeCharges:
    """
    Class to compute Ionicity and statistics from CHARGE.lobster data.

    :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
    :param path_to_charge: path to CHARGE.lobster (e.g., `CHARGE.lobster`)
    :param charge_type: set charge type used for computing ionicity.
        Possible options are `Mulliken` or `Loewdin`
    """

    def __init__(
        self,
        path_to_structure: str | Path,
        path_to_charge: str | Path,
        charge_type: Literal["mulliken", "loewdin"],
    ):
        """
        Compute the Ionicity of the structure and charge statistics from CHARGE.lobster data.

        :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
        :param path_to_charge: path to CHARGE.lobster (e.g., `CHARGE.lobster`)
        :param charge_type: set charge type used for computing ionicity.
            Possible options are `mulliken` or `loewdin`
        """
        self.path_to_structure = path_to_structure
        self.path_to_charge = path_to_charge
        self.charge_type = charge_type

        if self.charge_type.lower() not in ["mulliken", "loewdin"]:
            raise ValueError("Please check the requested charge_type. Possible options are `mulliken` or `loewdin`")

    def _calc_ionicity(self) -> float:
        r"""
        Calculate the ionicity of the crystal structure based on quantum chemical charges.

        References:
            - R. Nelson, C. Ertural, P. C. Müller, R. Dronskowski, 2023, DOI 10.1016/B978-0-12-823144-9.00120-5

        Returns:
            Ionicity of the structure

        """
        chargeobj = Charge(filename=self.path_to_charge)
        structure = Structure.from_file(self.path_to_structure)

        ch_veff = []
        tol = 1e-6
        for i, j in enumerate(getattr(chargeobj, self.charge_type.capitalize())):
            if (
                j > tol
                and not structure.species[i].is_transition_metal
                and (not structure.species[i].is_actinoid and not structure.species[i].is_lanthanoid)
            ):
                valence_elec = element(structure.species[i].value)
                try:
                    val = j / valence_elec.nvalence()
                except ZeroDivisionError:
                    val = 0
                ch_veff.append(val)

            elif (
                j < tol
                and not structure.species[i].is_transition_metal
                and (not structure.species[i].is_actinoid and not structure.species[i].is_lanthanoid)
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 8) if valence_elec.nvalence() != 8 else 0
                ch_veff.append(val)

            elif j > tol and (
                structure.species[i].is_transition_metal
                or structure.species[i].is_actinoid
                or structure.species[i].is_lanthanoid
            ):
                val = (
                    j / (structure.species[i].max_oxidation_state - 0)
                    if structure.species[i].max_oxidation_state != tol
                    else 0
                )
                ch_veff.append(val)

            elif j < tol and (
                structure.species[i].is_transition_metal
                or structure.species[i].is_actinoid
                or structure.species[i].is_lanthanoid
            ):
                val = (
                    j / (structure.species[i].max_oxidation_state - 8)
                    if structure.species[i].max_oxidation_state != 8
                    else 0
                )
                ch_veff.append(val)

        return sum(ch_veff) / structure.num_sites

    def _calc_stats(self) -> dict[str, float]:
        """
        Calculate standard statistics of the atomic-charges in CHARGE.lobster.

        Returns:
            A dictionary with charge statistics
        """
        chargeobj = Charge(filename=self.path_to_charge)
        charges = getattr(chargeobj, self.charge_type.capitalize())
        return {
            f"{self.charge_type.capitalize()}_mean": np.mean(charges),
            f"{self.charge_type.capitalize()}_min": np.min(charges),
            f"{self.charge_type.capitalize()}_max": np.max(charges),
            f"{self.charge_type.capitalize()}_std": np.std(charges),
        }

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        Return a pandas dataframe with computed ionicity and charge statistics as columns.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.

        Returns:
            Returns a pandas dataframe with ionicity and charge statistics as columns.

        """
        if not ids:
            ids = Path(self.path_to_charge).parent.name

        data = self._calc_stats()

        if self.charge_type.lower() == "mulliken":
            data["Ionicity_Mull"] = self._calc_ionicity()
        else:
            data["Ionicity_Loew"] = self._calc_ionicity()

        return pd.DataFrame(index=[ids], data=data)


class FeaturizeDoscar:
    """
    Class to compute DOS moments and fingerprints from DOSCAR.lobster / DOSCAR.LSO.lobster.

    :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
    :param path_to_doscar: path to DOSCAR.lobster or DOSCAR.LSO.lobster
    :param e_range: range of energy relative to fermi for which moment features and
        features needs to be computed
    :param add_element_dos_moments: add element dos moment features alongside orbital dos
    """

    def __init__(
        self,
        path_to_structure: str | Path,
        path_to_doscar: str | Path,
        e_range: list[float] | None = [-10.0, 0.0],
        add_element_dos_moments: bool = False,
    ):
        """
        Featurize DOSCAR.lobster or DOSCAR.LSO.lobster data.

        :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
        :param path_to_doscar: path to DOSCAR.lobster or DOSCAR.LSO.lobster
        :param e_range: range of energy relative to fermi for which moment features and
            features needs to be computed
        :param add_element_dos_moments: add element dos moment features alongside orbital dos
        """
        self.path_to_structure = path_to_structure
        self.path_to_doscar = path_to_doscar
        self.e_range = e_range
        self.dos = Doscar(doscar=self.path_to_doscar, structure_file=self.path_to_structure).completedos
        self.add_element_dos_moments = add_element_dos_moments

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        Return a pandas dataframe with computed DOS moment features as columns.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.

        Moment features are PDOS center, width, skewness, kurtosis and upper band edge.

        Returns:
            A pandas dataframe object
        """
        if ids:
            df = pd.DataFrame(index=[ids])
        else:
            ids = Path(self.path_to_doscar).parent.name
            df = pd.DataFrame(index=[ids])

        spd_dos_lobster = self.dos.get_spd_dos()

        for orbital in spd_dos_lobster:
            df.loc[ids, f"{orbital.name}_band_center"] = round(
                self.dos.get_band_center(band=orbital, erange=self.e_range), 4
            )
            df.loc[ids, f"{orbital.name}_band_width"] = round(
                self.dos.get_band_width(band=orbital, erange=self.e_range), 4
            )
            df.loc[ids, f"{orbital.name}_band_skew"] = round(
                self.dos.get_band_skewness(band=orbital, erange=self.e_range), 4
            )
            df.loc[ids, f"{orbital.name}_band_kurtosis"] = round(
                self.dos.get_band_kurtosis(band=orbital, erange=self.e_range), 4
            )
            df.loc[ids, f"{orbital.name}_band_upperband_edge"] = round(
                self.dos.get_upper_band_edge(band=orbital, erange=self.e_range), 4
            )
            if self.add_element_dos_moments:
                elements_list = self.dos.structure.elements
                for el in elements_list:
                    if orbital in self.dos.get_element_spd_dos(el=el.name):
                        df.loc[ids, f"{el.name}_{orbital.name}_band_center"] = round(
                            self.dos.get_band_center(band=orbital, elements=[el], erange=self.e_range),
                            4,
                        )
                        df.loc[ids, f"{el.name}_{orbital.name}_band_width"] = round(
                            self.dos.get_band_width(band=orbital, elements=[el], erange=self.e_range),
                            4,
                        )
                        df.loc[ids, f"{el.name}_{orbital.name}_band_skew"] = round(
                            self.dos.get_band_skewness(band=orbital, elements=[el], erange=self.e_range),
                            4,
                        )
                        df.loc[ids, f"{el.name}_{orbital.name}_band_kurtosis"] = round(
                            self.dos.get_band_kurtosis(band=orbital, elements=[el], erange=self.e_range),
                            4,
                        )
                        df.loc[ids, f"{el.name}_{orbital.name}_band_upperband_edge"] = round(
                            self.dos.get_upper_band_edge(band=orbital, elements=[el], erange=self.e_range),
                            4,
                        )

        return df

    def get_fingerprint_df(
        self,
        ids: str | None = None,
        fp_type: Literal["s", "p", "d", "f", "summed_pdos", "tdos"] = "summed_pdos",
        binning: bool = True,
        n_bins: int = 256,
        normalize: bool = True,
    ) -> pd.DataFrame:
        """
        Generate a dataframe consisting of DOS fingerprint (fp).

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.
        :param fp_type: Specify fingerprint type to compute, can accept `s/p/d/f/tdos/summed_pdos`
            (default is summed_pdos)
        :param binning: If true, the DOS fingerprint is binned using np.linspace and n_bins.
            Default is True.
        :param n_bins: Number of bins to be used in the fingerprint (default is 256)
        :param normalize: If true, normalizes the area under fp to equal to 1. Default is True.

        Returns:
            A pandas dataframe object with DOS fingerprints
        """
        if ids:
            df = pd.DataFrame(index=[ids], columns=["DOS_FP"])
        else:
            ids = Path(self.path_to_doscar).parent.name
            df = pd.DataFrame(index=[ids], columns=["DOS_FP"])

        dos_fp = self.dos.get_dos_fp(
            fp_type=fp_type,
            normalize=normalize,
            n_bins=n_bins,
            binning=binning,
            max_e=self.e_range[-1] if self.e_range is not None else None,
            min_e=self.e_range[0] if self.e_range is not None else None,
        )

        df.loc[ids, "DOS_FP"] = dos_fp

        return df


class FeaturizeIcoxxlist:
    """
    Class to Featurize ICOXXLIST.lobster as Bond weighted distribution function (BWDF).

    :param path_to_icoxxlist: path to ICOXXLIST.lobster
    :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
    :param path_to_grosspop: path to GROSSPOP.lobster
    :param bin_width: bin width for the BWDF
    :param interactions_tol: tolerance for interactions
    :param max_length: maximum bond length for BWDF computation
    :param min_length: minimum bond length for BWDF computation
    :param n_electrons_scaling: bool indicating if ICOXX values should be scaled by number of electrons.
        Only for testing purposes. Should not affect the results in any meaningful way.
    :param normalization: normalization strategy for BWDF
    :param are_cobis: bool indicating if file contains COBI/ICOBI data
    :param are_coops: bool indicating if file contains COOP/ICOOP data
    """

    def __init__(
        self,
        path_to_icoxxlist: str | Path,
        path_to_structure: str | Path,
        path_to_grosspop: str | Path | None = None,
        bin_width: float = 0.02,
        interactions_tol: float = 1e-3,
        max_length: float = 6.0,
        min_length: float = 0.0,
        n_electrons_scaling: bool = False,
        normalization: Literal["formula_units", "area", "counts", "none"] = "formula_units",
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
        """
        Initialize FeaturizeIcoxxlist class attributes.

        :param path_to_icoxxlist: path to ICOXXLIST.lobster
        :param path_to_structure: path to structure file (e.g., `CONTCAR` (preferred), `POSCAR`)
        :param path_to_grosspop: path to GROSSPOP.lobster
        :param bin_width: bin width for the BWDF
        :param interactions_tol: numerical tolerance considered for interactions to be insignificant
        :param max_length: maximum bond length for BWDF computation
        :param min_length: minimum bond length for BWDF computation
        :param n_electrons_scaling: bool indicating if ICOXX values should be scaled by number of electrons
        :param normalization: normalization strategy for BWDF
        :param are_cobis: bool indicating if file contains COBI/ICOBI data
        :param are_coops: bool indicating if file contains COOP/ICOOP data
        """
        self.path_to_icoxxlist = path_to_icoxxlist
        self.path_to_structure = path_to_structure
        self.path_to_grosspop = path_to_grosspop
        self.bin_width = bin_width
        self.are_cobis = are_cobis
        self.are_coops = are_coops
        self.icoxxlist = Icohplist(
            filename=self.path_to_icoxxlist,
            are_cobis=self.are_cobis,
            are_coops=self.are_coops,
        )
        self.grosspop = Grosspop(filename=self.path_to_grosspop) if self.path_to_grosspop else None
        self.interactions_tol = interactions_tol
        self.structure = Structure.from_file(self.path_to_structure)
        self.max_length = max_length
        self.min_length = min_length
        self.n_electrons_scaling = n_electrons_scaling
        self.normalization = normalization

    def get_icoxx_neighbors_data(self, site_index: int | None = None) -> dict:
        """
        Get the neighbors data with icoxx values for a structure.

        Uses a distance based neighbor list as reference to map the neighbor's data.

        Args:
            site_index: index of the site for which neighbors data is returned. Default is None (All sites).

        Returns:
            dict:
            Neighbors data as a dictionary with the following information

            - "ref_rdf_data": radial distribution function (RDF) data
            - "input_icoxx_list": complete ICOXXLIST.lobster data in the form of list of tuples
            - "mapped_icoxx_data": ICOXX values mapped to RDF data
            - "missing_interactions": list of interactions that are present in RDF data but not in ICOXX data
            - "wasserstein_dist_to_rdf": wasserstein distance computed between ref_rdf_data and mapped_icoxx_data.
        """
        # Get all neighbors data in a single list using distance based algorithm
        # This is used as a reference to map the icoxx data
        rdf_nb_lst = self.structure.get_neighbor_list(r=self.max_length)

        # Collect bond lengths, atom labels and bond strengths
        if site_index is not None:
            if site_index not in range(self.structure.num_sites):
                raise ValueError(f"{site_index} is not a valid site index for the structure")
            bond_labels = np.array(self.icoxxlist.icohpcollection._list_labels)
            indices = np.where(
                np.isin(
                    bond_labels,
                    list(
                        self.icoxxlist.icohpcollection.get_icohp_dict_of_site(
                            site_index, minbondlength=self.min_length, maxbondlength=self.max_length
                        ).keys()
                    ),
                )
            )[0]

            bond_lengths = np.array(self.icoxxlist.icohpcollection._list_length)[indices].tolist()
            atoms1 = np.array(self.icoxxlist.icohpcollection._list_atom1)[indices].tolist()
            atoms2 = np.array(self.icoxxlist.icohpcollection._list_atom2)[indices].tolist()
            icoxxs = np.array([sum(item.values()) for item in self.icoxxlist.icohpcollection._list_icohp])[
                indices
            ].tolist()
            trans = [list(item) for item in np.array(self.icoxxlist.icohpcollection._list_translation)[indices]]
            pairs = [[at1, at2] for at1, at2 in zip(atoms1, atoms2)]

            if self.n_electrons_scaling and self.grosspop:
                atom_1_index = np.array([int("".join(filter(str.isdigit, at1))) - 1 for at1 in atoms1])
                atom_2_index = np.array([int("".join(filter(str.isdigit, at2))) - 1 for at2 in atoms2])
                pair_n_electrons = np.array(
                    [
                        self.grosspop.list_dict_grosspop[at1]["Loewdin GP"]["total"]
                        + self.grosspop.list_dict_grosspop[at2]["Loewdin GP"]["total"]
                        for at1, at2 in zip(atom_1_index, atom_2_index)
                    ]
                )
                icoxxs = (np.array(icoxxs) / pair_n_electrons).tolist()

            relevant_site_indices = [
                ix for ix, (org, dest) in enumerate(zip(rdf_nb_lst[0], rdf_nb_lst[1])) if org == site_index
            ]
            rdf_distances = [round(dist, 5) for ix, dist in enumerate(rdf_nb_lst[3]) if ix in relevant_site_indices]
            rdf_trans = [[int(i) for i in img] for ix, img in enumerate(rdf_nb_lst[2]) if ix in relevant_site_indices]
            rdf_pairs = [
                [self.structure[org].species_string + str(org + 1), self.structure[dest].species_string + str(dest + 1)]
                for ix, (org, dest) in enumerate(zip(rdf_nb_lst[0], rdf_nb_lst[1]))
                if ix in relevant_site_indices
            ]

        else:  # complete structure
            bond_lengths = self.icoxxlist.icohpcollection._list_length
            atoms1 = self.icoxxlist.icohpcollection._list_atom1
            atoms2 = self.icoxxlist.icohpcollection._list_atom2
            icoxxs = [sum(item.values()) for item in self.icoxxlist.icohpcollection._list_icohp]
            trans = [list(item) for item in self.icoxxlist.icohpcollection._list_translation]
            pairs = [[at1, at2] for at1, at2 in zip(atoms1, atoms2)]

            if self.n_electrons_scaling and self.grosspop:
                atom_1_index = np.array([int("".join(filter(str.isdigit, at1))) - 1 for at1 in atoms1])
                atom_2_index = np.array([int("".join(filter(str.isdigit, at2))) - 1 for at2 in atoms2])
                pair_n_electrons = np.array(
                    [
                        self.grosspop.list_dict_grosspop[at1]["Loewdin GP"]["total"]
                        + self.grosspop.list_dict_grosspop[at2]["Loewdin GP"]["total"]
                        for at1, at2 in zip(atom_1_index, atom_2_index)
                    ]
                )
                icoxxs = (np.array(icoxxs) / pair_n_electrons).tolist()

            rdf_distances = [round(dist, 5) for dist in rdf_nb_lst[3]]
            rdf_trans = [[int(i) for i in img] for img in rdf_nb_lst[2]]
            rdf_pairs = [
                [self.structure[org].species_string + str(org + 1), self.structure[dest].species_string + str(dest + 1)]
                for org, dest in zip(rdf_nb_lst[0], rdf_nb_lst[1])
            ]

        # Assimilate icoxx data in tuples
        input_icoxx_list = list(zip(pairs, bond_lengths, trans))

        # 0: pair, 1: bond length, 2: translation
        ref_rdf_data = list(zip(rdf_pairs, rdf_distances, rdf_trans))

        # Create a temp dict for faster lookup
        icoxx_dict = {
            (tuple(icoxx_p_d_t[0]), tuple(icoxx_p_d_t[2])): (icoxx_p_d_t[1], icoxxs[idx])
            for idx, icoxx_p_d_t in enumerate(input_icoxx_list)
        }
        # Initialize an empty list to store the data that goes in BWDF computation
        mapped_icoxx_data, missing_interactions = [], []

        # Check if interaction exists in both rdf data / icoxx data, then add to complete data
        for rdf_p_d_t in ref_rdf_data:
            pair_translation_key = (tuple(rdf_p_d_t[0]), tuple(rdf_p_d_t[2]))
            bond_length = rdf_p_d_t[1]
            if pair_translation_key in icoxx_dict:
                icoxx_bond_length, icoxx_value = icoxx_dict[pair_translation_key]
                if round(np.absolute(bond_length - icoxx_bond_length), 5) <= 1e-5:
                    complete_entry = (*rdf_p_d_t, icoxx_value)
                    if complete_entry not in mapped_icoxx_data:
                        mapped_icoxx_data.append(complete_entry)
                else:
                    missing_interactions.append(rdf_p_d_t)
            else:
                missing_interactions.append(rdf_p_d_t)

        # Check if the missing interactions are reverse images in ref neighbours data which LOBSTER sometimes eliminates
        for rdf_p_d_t in missing_interactions[:]:
            reversed_translation = tuple(-1 * np.array(rdf_p_d_t[2]))
            reversed_pair = tuple(reversed(rdf_p_d_t[0]))
            reversed_key = (reversed_pair, reversed_translation)

            if reversed_key in icoxx_dict:
                icoxx_bond_length, icoxx_value = icoxx_dict[reversed_key]
                if round(np.absolute(rdf_p_d_t[1] - icoxx_bond_length), 5) <= 1e-5:
                    complete_entry = (*rdf_p_d_t, icoxx_value)
                    if complete_entry not in mapped_icoxx_data:
                        mapped_icoxx_data.append((*rdf_p_d_t, icoxx_value))
                        missing_interactions.remove(rdf_p_d_t)

        # Compute histogram of bond lengths from the reference rdf and
        # icoxx data to get wasserstein distance
        rdf_dist = [entry[1] for entry in ref_rdf_data]
        icoxx_dist = [entry[1] for entry in mapped_icoxx_data]
        rdf_hist, _ = np.histogram(
            rdf_dist,
            bins=np.arange(0, self.max_length + self.bin_width, self.bin_width),
            density=False,
        )
        icoxx_hist, _ = np.histogram(
            icoxx_dist,
            bins=np.arange(0, self.max_length + self.bin_width, self.bin_width),
            density=False,
        )

        return {
            "ref_rdf_data": ref_rdf_data,
            "input_icoxx_list": input_icoxx_list,
            "mapped_icoxx_data": mapped_icoxx_data,
            "missing_interactions": missing_interactions,
            "wasserstein_dist_to_rdf": wasserstein_distance(rdf_hist, icoxx_hist),
        }

    def calc_bwdf(self) -> dict:
        """
        Compute BWDF from ICOXXLIST.lobster data.

        Returns:
            dict:
            BWDF as a dictionary for each atom pair and entire structure

            - "A-B": BWDF for atom pair A-B, e.g., "Na-Cl": {“icoxx_binned”: np.array, “icoxx_counts”: np.array}
            - "summed": BWDF for entire structure, e.g., "summed": {“icoxx_binned”: np.array, “icoxx_counts”: np.array}
            - "centers": bin centers for BWDF
            - "edges": bin edges for BWDF
            - "bin_width": bin width
            - "wasserstein_dist_to_rdf": wasserstein distance between RDF and ICOXX data

        """
        # Get all neighbors data in a single list
        icoxx_neighbors_data = self.get_icoxx_neighbors_data()

        # Extract unique atomic species and create possible combinations
        species_combinations = [sorted(item) for item in combinations_with_replacement(self.structure.symbol_set, 2)]

        # Calculate number of bins
        n_bins = int(np.ceil(((self.max_length + self.bin_width) - self.min_length) / self.bin_width))

        # Get bin edges and centers
        bin_edges = np.round(np.linspace(self.min_length, self.max_length, n_bins), 5)
        bin_centers = bin_edges[:-1] + self.bin_width / 2

        # Initialize dictionary for storing binned data by atom pair
        bwdf_atom_pair = {
            "-".join(atom_pair): {
                "icoxx_binned": np.zeros(bin_centers.shape),
                "icoxx_counts": np.zeros(bin_centers.shape),
            }
            for atom_pair in species_combinations
        }

        # Populate bins with corresponding icoxx values
        for key in bwdf_atom_pair:
            for interactions in icoxx_neighbors_data["mapped_icoxx_data"]:
                sorted_entry = sorted([interactions[0][0].strip("0123456789"), interactions[0][1].strip("0123456789")])
                if key == "-".join(sorted_entry):
                    for ii, l1, l2 in zip(range(len(bin_centers)), bin_edges[:-1], bin_edges[1:], strict=False):
                        if l1 <= interactions[1] < l2:
                            # Add icoxx values to the corresponding bin
                            bwdf_atom_pair[key]["icoxx_binned"][ii] += (
                                interactions[3] if abs(interactions[3]) > self.interactions_tol else 0
                            )
                            bwdf_atom_pair[key]["icoxx_counts"][ii] += (
                                1 if abs(interactions[3]) > self.interactions_tol else 0
                            )

        icoxx_binned_summed = np.sum(
            [bwdf_atom_pair[atom_pair]["icoxx_binned"] for atom_pair in bwdf_atom_pair], axis=0
        )
        icoxx_counts_summed = np.sum(
            [bwdf_atom_pair[atom_pair]["icoxx_counts"] for atom_pair in bwdf_atom_pair], axis=0
        )

        bwdf_atom_pair["summed"] = {"icoxx_binned": icoxx_binned_summed, "icoxx_counts": icoxx_counts_summed}
        bwdf_atom_pair["centers"] = bin_centers
        bwdf_atom_pair["edges"] = bin_edges
        bwdf_atom_pair["bin_width"] = self.bin_width  # type: ignore[assignment]
        bwdf_atom_pair["wasserstein_dist_to_rdf"] = icoxx_neighbors_data["wasserstein_dist_to_rdf"]

        # Normalize BWDF data
        return self._normalize_bwdf(bwdf=bwdf_atom_pair)

    def calc_site_asymmetry_index(self, site_index: int) -> float:
        """
        Compute the asymmetry index for a site using bond strengths as weights.

        Args:
            site_index: index of the site for which the asymmetry index needs to be computed

        References:
            - F. Belli, E. Zurek, I. Errea, 2025, DOI 10.48550/arXiv.2501.14420

        Returns:
            Asymmetry index for the site
        """
        if site_index not in range(self.structure.num_sites):
            raise ValueError(f"{site_index} is not a valid site index for the structure")

        # Get all neighbors data in a single list
        icoxx_neighbors_data = self.get_icoxx_neighbors_data(site_index=site_index)["mapped_icoxx_data"]

        # Calculate the asymmetry index
        icoxxs_x_y_z = []
        for pair in icoxx_neighbors_data:
            src_dst = [int("".join(filter(str.isdigit, p))) - 1 for p in pair[0]]
            src = src_dst[0]
            dst = src_dst[1]
            cart_dst = self.structure.lattice.get_cartesian_coords(self.structure[dst].frac_coords + np.array(pair[2]))
            cart_src = self.structure[src].coords
            unit_vec = (cart_dst - cart_src) / pair[1]
            icoxxs_x_y_z.append(pair[-1] * unit_vec)

        return np.linalg.norm(np.mean(icoxxs_x_y_z, axis=0))

    def calc_site_bwdf(self, site_index: int) -> dict:
        """
        Compute BWDF from ICOXXLIST.lobster data for a site.

        Args:
            site_index: index of the site for which BWDF needs to be computed

        Returns:
            dict:
            BWDF as a dictionary for the site in the following format

            - "X": BWDF for the site X, e.g., "0": {“icoxx_binned”: np.array, “icoxx_counts”: np.array}
            - "centers": bin centers for BWDF
            - "edges": bin edges for BWDF
            - "bin_width": bin width
            - "wasserstein_dist_to_rdf": wasserstein distance between RDF and ICOXX data

        """
        if site_index not in range(self.structure.num_sites):
            raise ValueError(f"{site_index} is not a valid site index for the structure")

        # Get all neighbors data in a single list
        icoxx_neighbors_data = self.get_icoxx_neighbors_data(site_index=site_index)

        # Calculate number of bins
        n_bins = int(np.ceil(((self.max_length + self.bin_width) - self.min_length) / self.bin_width))

        # Get bin edges and centers
        bin_edges = np.round(np.linspace(self.min_length, self.max_length, n_bins), 5)
        bin_centers = bin_edges[:-1] + self.bin_width / 2

        # Initialize dictionary for storing binned data by atom pair
        site_bwdf = {
            f"{site_index}": {"icoxx_binned": np.zeros(bin_centers.shape), "icoxx_counts": np.zeros(bin_centers.shape)}
        }

        for interactions in icoxx_neighbors_data["mapped_icoxx_data"]:
            for ii, l1, l2 in zip(range(len(bin_centers)), bin_edges[:-1], bin_edges[1:]):
                if l1 <= interactions[1] < l2:
                    # Add icoxx values to the corresponding bin
                    site_bwdf[f"{site_index}"]["icoxx_binned"][ii] += (
                        interactions[3] if abs(interactions[3]) > self.interactions_tol else 0
                    )
                    site_bwdf[f"{site_index}"]["icoxx_counts"][ii] += (
                        1 if abs(interactions[3]) > self.interactions_tol else 0
                    )

        site_bwdf["centers"] = bin_centers
        site_bwdf["edges"] = bin_edges
        site_bwdf["bin_width"] = self.bin_width  # type: ignore[assignment]
        site_bwdf["wasserstein_dist_to_rdf"] = icoxx_neighbors_data["wasserstein_dist_to_rdf"]  # type: ignore[assignment]

        # Normalize BWDF data
        return self._normalize_bwdf(bwdf=site_bwdf)

    def calc_label_bwdf(self, bond_label: str) -> dict:
        """
        Compute BWDF from ICOXXLIST.lobster data for a bond label.

        Args:
            bond_label: bond label for which BWDF needs to be computed

        Returns:
            dict:
            BWDF as a dictionary for the bond label in the following format

            - "X": BWDF for the bond label, e.g., "20": {“icoxx_binned”: np.array, “icoxx_counts”: np.array}
            - "centers": bin centers for BWDF
            - "edges": bin edges for BWDF
            - "bin_width": bin width
            - "wasserstein_dist_to_rdf": wasserstein distance between RDF and ICOXX data
        """
        index = self.icoxxlist.icohpcollection._list_labels.index(bond_label)
        bond_length = self.icoxxlist.icohpcollection._list_length[index]
        atom1 = self.icoxxlist.icohpcollection._list_atom1[index]
        atom2 = self.icoxxlist.icohpcollection._list_atom2[index]
        icoxx = sum(self.icoxxlist.icohpcollection._list_icohp[index].values())
        trans = self.icoxxlist.icohpcollection._list_translation[index]

        # Complete data
        complete_data = [(sorted([atom1, atom2]), bond_length, trans, icoxx)]

        # Calculate number of bins
        n_bins = int(np.ceil(((self.max_length + self.bin_width) - self.min_length) / self.bin_width))

        # Get bin edges and centers
        bin_edges = np.round(np.linspace(self.min_length, self.max_length, n_bins), 5)
        bin_centers = bin_edges[:-1] + self.bin_width / 2

        # Initialize dictionary for storing binned data by atom pair
        label_bwdf = {
            bond_label: {"icoxx_binned": np.zeros(bin_centers.shape), "icoxx_counts": np.zeros(bin_centers.shape)}
        }

        for interactions in complete_data:
            for ii, l1, l2 in zip(range(len(bin_centers)), bin_edges[:-1], bin_edges[1:]):
                if l1 <= interactions[1] < l2:
                    # Add icoxx values to the corresponding bin
                    label_bwdf[bond_label]["icoxx_binned"][ii] += (
                        interactions[3] if abs(interactions[3]) > self.interactions_tol else 0
                    )
                    label_bwdf[bond_label]["icoxx_counts"][ii] += (
                        1 if abs(interactions[3]) > self.interactions_tol else 0
                    )
        label_bwdf["centers"] = bin_centers
        label_bwdf["edges"] = bin_edges
        label_bwdf["bin_width"] = self.bin_width  # type: ignore[assignment]

        # Normalize BWDF data
        return self._normalize_bwdf(bwdf=label_bwdf)

    @staticmethod
    def _get_features_col_names(bwdf: dict) -> list[str]:
        """
        Get the names of the features that will be generated by the class.

        Args:
            bwdf: dictionary containing the BWDF data

        Returns:
            list of feature names
        """
        features = []
        for edge_1, edge_2 in zip(bwdf["edges"][:-1], bwdf["edges"][1:]):
            features.append(f"bwdf_{round(edge_1, 2)}-{round(edge_2, 2)}")

        return features

    def _normalize_bwdf(self, bwdf: dict) -> dict:
        """
        Normalize BWDF data based on the normalization strategy.

        :param bwdf: BWDF data as a dictionary

        Returns:
            Normalized BWDF data as a dictionary
        """
        for bwdf_label, bwdf_value in bwdf.items():
            if bwdf_label not in ("centers", "edges", "bin_width", "wasserstein_dist_to_rdf"):
                if self.normalization == "area":
                    total_area = np.sum(np.abs(bwdf_value["icoxx_binned"]) * self.bin_width)
                    bwdf[bwdf_label]["icoxx_binned"] = np.nan_to_num(bwdf_value["icoxx_binned"] / total_area)
                elif self.normalization == "formula_units":
                    formula_units = self.structure.composition.get_reduced_formula_and_factor()[-1]
                    bwdf[bwdf_label]["icoxx_binned"] = np.nan_to_num(bwdf_value["icoxx_binned"] / formula_units)
                elif self.normalization == "counts":
                    bwdf[bwdf_label]["icoxx_binned"] = np.nan_to_num(
                        bwdf_value["icoxx_binned"] / bwdf_value["icoxx_counts"]
                    )

        return bwdf

    def get_asymmetry_index_stats_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with asymmetry index statistical information as columns.

        Args:
              ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with asymmetry index statistical information as columns.
            Columns include sum, mean, std, min, and max.
        """
        asymmetry_indices = []

        for site_index in range(self.structure.num_sites):
            asymmetry_indices.append(self.calc_site_asymmetry_index(site_index))

        column_names = ["asi_sum", "asi_mean", "asi_std", "asi_min", "asi_max"]

        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        df.loc[ids, "asi_sum"] = np.sum(asymmetry_indices)
        df.loc[ids, "asi_mean"] = np.mean(asymmetry_indices)
        df.loc[ids, "asi_std"] = np.std(asymmetry_indices)
        df.loc[ids, "asi_min"] = np.min(asymmetry_indices)
        df.loc[ids, "asi_max"] = np.max(asymmetry_indices)

        return df

    def get_binned_bwdf_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with computed BWDF features as columns.

        Args:
            ids: set index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with BWDF as columns. Each column contains
            sum of icoxx values corresponding to bins.
        """
        bwdf = self.calc_bwdf()
        column_names = self._get_features_col_names(bwdf=bwdf)
        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        for icoxx_weight, col in zip(bwdf["summed"]["icoxx_binned"], column_names):
            df.loc[ids, col] = icoxx_weight

        return df

    def get_site_df(self, site_index: int, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with computed BWDF features for a site as columns.

        Args:
            site_index: index of the site in a structure for which BWDF needs to be computed
            ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with BWDF as columns for a site. Each column contains
            sum of icoxx values corresponding to bins.
        """
        site_bwdf = self.calc_site_bwdf(site_index=site_index)
        column_names = [f"{feat_name}_site_{site_index}" for feat_name in self._get_features_col_names(bwdf=site_bwdf)]
        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        for icoxx_weight, col in zip(site_bwdf[f"{site_index}"]["icoxx_binned"], column_names):
            df.loc[ids, col] = icoxx_weight

        return df

    def get_site_bwdf_stats_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas datafram with mean and std from sitewise BWDFs.

        Args:
            ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with BWDF statistical information as columns.
            The columns include the mean and standard deviation calculated from the
            sitewise BWDFs stats (i.e., sum, mean, minimum, maximum, std, skewness, and kurtosis).
        """
        column_names = [
            "site_bwdf_sum_mean",
            "site_bwdf_mean_mean",
            "site_bwdf_std_mean",
            "site_bwdf_min_mean",
            "site_bwdf_max_mean",
            "site_bwdf_skew_mean",
            "site_bwdf_kurtosis_mean",
            "site_bwdf_sum_std",
            "site_bwdf_mean_std",
            "site_bwdf_std_std",
            "site_bwdf_min_std",
            "site_bwdf_max_std",
            "site_bwdf_skew_std",
            "site_bwdf_kurtosis_std",
        ]
        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        bwdf_sums = []
        bwdf_means = []
        bwdf_stds = []
        bwdf_mins = []
        bwdf_maxs = []
        bwdf_skews = []
        bwdf_kurtosis = []
        for site_index in range(self.structure.num_sites):
            site_bwdf = self.calc_site_bwdf(site_index=site_index)

            bwdf_sums.append(
                np.sum(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(np.sum(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_means.append(
                np.mean(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(np.mean(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_stds.append(
                np.std(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(np.std(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_mins.append(
                np.min(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(np.min(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_maxs.append(
                np.max(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(np.max(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_skews.append(
                skew(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(skew(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )
            bwdf_kurtosis.append(
                kurtosis(site_bwdf[f"{site_index}"]["icoxx_binned"])
                if not np.isnan(kurtosis(site_bwdf[f"{site_index}"]["icoxx_binned"]))
                else 0
            )

        df.loc[ids, "site_bwdf_sum_mean"] = np.mean(bwdf_sums)
        df.loc[ids, "site_bwdf_mean_mean"] = np.mean(bwdf_means)
        df.loc[ids, "site_bwdf_std_mean"] = np.mean(bwdf_stds)
        df.loc[ids, "site_bwdf_min_mean"] = np.mean(bwdf_mins)
        df.loc[ids, "site_bwdf_max_mean"] = np.mean(bwdf_maxs)
        df.loc[ids, "site_bwdf_skew_mean"] = np.mean(bwdf_skews)
        df.loc[ids, "site_bwdf_kurtosis_mean"] = np.mean(bwdf_kurtosis)
        df.loc[ids, "site_bwdf_sum_std"] = np.std(bwdf_sums)
        df.loc[ids, "site_bwdf_mean_std"] = np.std(bwdf_means)
        df.loc[ids, "site_bwdf_std_std"] = np.std(bwdf_stds)
        df.loc[ids, "site_bwdf_min_std"] = np.std(bwdf_mins)
        df.loc[ids, "site_bwdf_max_std"] = np.std(bwdf_maxs)
        df.loc[ids, "site_bwdf_skew_std"] = np.std(bwdf_skews)
        df.loc[ids, "site_bwdf_kurtosis_std"] = np.std(bwdf_kurtosis)

        return df

    def get_pair_bwdf_stats_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with statistical info from pairwise BWDFs.

        Args:
            ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with BWDF statistical information as columns.
            The columns include the mean and standard deviation calculated from the
            pairwise BWDFs stats (i.e., sum, mean, minimum, maximum, std, skewness, and kurtosis).
        """
        column_names = [
            "pair_bwdf_sum_mean",
            "pair_bwdf_mean_mean",
            "pair_bwdf_std_mean",
            "pair_bwdf_min_mean",
            "pair_bwdf_max_mean",
            "pair_bwdf_skew_mean",
            "pair_bwdf_kurtosis_mean",
            "pair_bwdf_sum_std",
            "pair_bwdf_mean_std",
            "pair_bwdf_std_std",
            "pair_bwdf_min_std",
            "pair_bwdf_max_std",
            "pair_bwdf_skew_std",
            "pair_bwdf_kurtosis_std",
        ]
        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        bwdf_sums = []
        bwdf_means = []
        bwdf_stds = []
        bwdf_mins = []
        bwdf_maxs = []
        bwdf_skews = []
        bwdf_kurtosis = []

        bwdf = self.calc_bwdf()
        for atom_pair in bwdf:
            if atom_pair not in ("summed", "centers", "edges", "bin_width", "wasserstein_dist_to_rdf"):
                bwdf_sums.append(
                    np.sum(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(np.sum(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )
                bwdf_means.append(
                    np.mean(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(np.mean(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )
                bwdf_stds.append(
                    np.std(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(np.std(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )
                bwdf_mins.append(
                    np.min(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(np.min(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )
                bwdf_maxs.append(
                    np.max(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(np.max(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )
                bwdf_skews.append(
                    skew(bwdf[atom_pair]["icoxx_binned"]) if not np.isnan(skew(bwdf[atom_pair]["icoxx_binned"])) else 0
                )
                bwdf_kurtosis.append(
                    kurtosis(bwdf[atom_pair]["icoxx_binned"])
                    if not np.isnan(kurtosis(bwdf[atom_pair]["icoxx_binned"]))
                    else 0
                )

        df.loc[ids, "pair_bwdf_sum_mean"] = np.mean(bwdf_sums)
        df.loc[ids, "pair_bwdf_mean_mean"] = np.mean(bwdf_means)
        df.loc[ids, "pair_bwdf_std_mean"] = np.mean(bwdf_stds)
        df.loc[ids, "pair_bwdf_min_mean"] = np.mean(bwdf_mins)
        df.loc[ids, "pair_bwdf_max_mean"] = np.mean(bwdf_maxs)
        df.loc[ids, "pair_bwdf_skew_mean"] = np.mean(bwdf_skews)
        df.loc[ids, "pair_bwdf_kurtosis_mean"] = np.mean(bwdf_kurtosis)
        df.loc[ids, "pair_bwdf_sum_std"] = np.std(bwdf_sums)
        df.loc[ids, "pair_bwdf_mean_std"] = np.std(bwdf_means)
        df.loc[ids, "pair_bwdf_std_std"] = np.std(bwdf_stds)
        df.loc[ids, "pair_bwdf_min_std"] = np.std(bwdf_mins)
        df.loc[ids, "pair_bwdf_max_std"] = np.std(bwdf_maxs)
        df.loc[ids, "pair_bwdf_skew_std"] = np.std(bwdf_skews)
        df.loc[ids, "pair_bwdf_kurtosis_std"] = np.std(bwdf_kurtosis)

        return df

    def get_summed_bwdf_stats_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with statistical info from BWDF as columns.

        Args:
              ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with BWDF statistical information as columns.
            Columns include sum, mean, std, min, max, skew, kurtosis, weighted mean and weighted std.
        """
        bwdf = self.calc_bwdf()
        bin_weights = np.abs(bwdf["summed"]["icoxx_binned"] / np.sum(bwdf["summed"]["icoxx_binned"]))
        column_names = [
            "bwdf_sum",
            "bwdf_mean",
            "bwdf_std",
            "bwdf_min",
            "bwdf_max",
            "bwdf_skew",
            "bwdf_kurtosis",
            "bwdf_w_mean",
            "bwdf_w_std",
        ]
        if ids:
            df = pd.DataFrame(index=[ids], columns=column_names)
        else:
            ids = Path(self.path_to_icoxxlist).parent.name
            df = pd.DataFrame(index=[ids], columns=column_names)

        w_bwdf_mean = np.average(bwdf["summed"]["icoxx_binned"], weights=bin_weights)
        w_bwdf_std = np.sqrt(np.average((bwdf["summed"]["icoxx_binned"] - w_bwdf_mean) ** 2, weights=bin_weights))

        df.loc[ids, "bwdf_sum"] = np.sum(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_mean"] = np.mean(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_std"] = np.std(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_w_mean"] = w_bwdf_mean
        df.loc[ids, "bwdf_w_std"] = w_bwdf_std
        df.loc[ids, "bwdf_min"] = np.min(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_max"] = np.max(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_skew"] = skew(bwdf["summed"]["icoxx_binned"])
        df.loc[ids, "bwdf_kurtosis"] = kurtosis(bwdf["summed"]["icoxx_binned"])

        return df

    def get_stats_df(
        self, ids: str | None = None, stats_type: Literal["atompair", "site", "summed", "all"] = "summed"
    ) -> pd.DataFrame:
        """Convenience method to get a pandas dataframe with statistical info from BWDF as columns.

        Args:
              ids: set the index name in the pandas dataframe. Default is None.
              stats_type: type of BWDF stats to be returned. Default is "summed".

                - "atompair": compute stats from unique atom pairs BWDFs.
                - "site": compute stats from site BWDFs.
                - "summed": compute stats from structure BWDFs.
                - "all": concatenated dataframe from `atompair`, `site` and `summed` options.

        Returns:
            A pandas dataframe object with BWDF statistical information as columns.
        """
        if stats_type == "atompair":
            df = self.get_pair_bwdf_stats_df(ids=ids)
        elif stats_type == "site":
            df = self.get_site_bwdf_stats_df(ids=ids)
        elif stats_type == "summed":
            df = self.get_summed_bwdf_stats_df(ids=ids)
        else:
            df = pd.concat(
                [
                    self.get_pair_bwdf_stats_df(ids=ids),
                    self.get_site_bwdf_stats_df(ids=ids),
                    self.get_summed_bwdf_stats_df(ids=ids),
                ],
                axis=1,
            )
        return df

    def get_sorted_bwdf_df(self, ids: str | None = None) -> pd.DataFrame:
        """Return a pandas dataframe with BWDF values sorted by distances, ascending.

        Args:
            ids: set the index name in the pandas dataframe. Default is None.

        Returns:
            A pandas dataframe object with binned BWDF values sorted by distance.
        """
        if not ids:
            ids = Path(self.path_to_icoxxlist).parent.name
        icoxx = self.calc_bwdf()["summed"]["icoxx_binned"]
        sorted_icoxx = icoxx[icoxx != 0.0]
        column_names = [f"bwdf_at_dist{d}" for d in range(len(sorted_icoxx))]
        return pd.DataFrame.from_dict({ids: dict(zip(column_names, sorted_icoxx))}).T

    def get_sorted_dist_df(
        self, ids: str | None = None, mode: Literal["positive", "negative"] = "negative"
    ) -> pd.DataFrame:
        """Return a pandas dataframe with distances sorted by BWDF values
        (either only positive or negative),  sorted descending by absolute values.

        Args:
            ids: set the index name in the pandas dataframe. Default is None
            mode: must be in ("positive", "negative"), defines whether BWDF values above or
                below zero are considered for distance featurization.

        Returns:
            A pandas dataframe object with binned distances sorted by BWDF values.
        """
        if mode not in ["positive", "negative"]:
            raise ValueError("param mode must be in ('positive', 'negative')")
        if not ids:
            ids = Path(self.path_to_icoxxlist).parent.name
        bwdf = self.calc_bwdf()
        sorted_dist_ids = bwdf["summed"]["icoxx_binned"].argsort()

        if mode == "negative":
            sorted_dists = [bwdf["centers"][i] for i in sorted_dist_ids if bwdf["summed"]["icoxx_binned"][i] < 0]
        else:
            sorted_dists = [bwdf["centers"][i] for i in sorted_dist_ids[::-1] if bwdf["summed"]["icoxx_binned"][i] > 0]

        column_names = [f"dist_at_{mode[:3]}_bwdf{d}" for d in range(len(sorted_dists))]
        return pd.DataFrame.from_dict({ids: dict(zip(column_names, sorted_dists))}).T
