# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""This module defines classes to featurize Lobster data ready to be used for ML studies."""

from __future__ import annotations

import gzip
import json
import warnings
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from mendeleev import element
from numpy import ndarray
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster import Charge, Doscar, Icohplist, MadelungEnergies
from scipy.integrate import trapezoid
from scipy.signal import hilbert

from lobsterpy.cohp.analyze import Analysis

from . import get_file_paths

warnings.filterwarnings("ignore")


class FeaturizeLobsterpy:
    """
    Class to featurize lobsterpy data.

    :param path_to_lobster_calc: path to parent directory containing lobster calc outputs
    :param path_to_json: path to lobster lightweight json
    :param bonds: "all" or "cation-anion" bonds
    :param orbital_resolved: bool indicating whether LobsterPy analysis is performed orbital wise
    """

    def __init__(
        self,
        path_to_lobster_calc: str | Path | None = None,
        path_to_json: str | Path | None = None,
        orbital_resolved: bool = False,
        bonds: str = "all",
    ):
        """Initialize featurizer."""
        self.path_to_json = path_to_json
        self.path_to_lobster_calc = path_to_lobster_calc
        self.orbital_resolved = orbital_resolved
        self.bonds = bonds

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
            if not ids:
                ids = Path(self.path_to_json).name.split(".")[0]

        elif self.path_to_lobster_calc and not self.path_to_json:
            # get lobsterpy condensed bonding analysis data using get_lobsterpy_cba_dict method
            data = FeaturizeLobsterpy.get_lobsterpy_cba_dict(
                path_to_lobster_calc=self.path_to_lobster_calc,
                bonds=self.bonds,
                orbital_resolved=self.orbital_resolved,
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
                        icohp_mean.append(float(bond_data["ICOHP_mean"]))
                        icohp_sum.append(float(bond_data["ICOHP_sum"]))
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
                                    icohp_mean_orb_bndg.append(bond_data["orbital_data"][orb_pair]["ICOHP_mean"])
                                    icohp_sum_orb_bndg.append(bond_data["orbital_data"][orb_pair]["ICOHP_sum"])
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
                                    icohp_mean_orb_antibndg.append(bond_data["orbital_data"][orb_pair]["ICOHP_mean"])
                                    icohp_sum_orb_antibndg.append(bond_data["orbital_data"][orb_pair]["ICOHP_sum"])
                                    antibond_orb.append(
                                        bond_data["orbital_data"][orb_pair]["orb_contribution_perc_antibonding"]
                                    )

        # add ICOHP stats data (mean, min, max, standard deviation) as columns to the dataframe

        df.loc[ids, "Icohp_mean_avg"] = 0 if len(icohp_mean) == 0 else np.mean(icohp_mean)
        df.loc[ids, "Icohp_mean_max"] = 0 if len(icohp_mean) == 0 else np.max(icohp_mean)
        df.loc[ids, "Icohp_mean_min"] = 0 if len(icohp_mean) == 0 else np.min(icohp_mean)
        df.loc[ids, "Icohp_mean_std"] = 0 if len(icohp_mean) == 0 else np.std(icohp_mean)

        df.loc[ids, "Icohp_sum_avg"] = 0 if len(icohp_sum) == 0 else np.mean(icohp_sum)
        df.loc[ids, "Icohp_sum_max"] = 0 if len(icohp_sum) == 0 else np.max(icohp_sum)
        df.loc[ids, "Icohp_sum_min"] = 0 if len(icohp_sum) == 0 else np.min(icohp_sum)
        df.loc[ids, "Icohp_sum_std"] = 0 if len(icohp_sum) == 0 else np.std(icohp_sum)

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
            df.loc[ids, "Icohp_bndg_orb_mean_avg"] = (
                0 if len(icohp_mean_orb_bndg) == 0 else np.mean(icohp_mean_orb_bndg)
            )
            df.loc[ids, "Icohp_bndg_orb_mean_max"] = 0 if len(icohp_mean_orb_bndg) == 0 else np.max(icohp_mean_orb_bndg)
            df.loc[ids, "Icohp_bndg_orb_mean_min"] = 0 if len(icohp_mean_orb_bndg) == 0 else np.min(icohp_mean_orb_bndg)
            df.loc[ids, "Icohp_bndg_orb_mean_std"] = 0 if len(icohp_mean_orb_bndg) == 0 else np.std(icohp_mean_orb_bndg)

            df.loc[ids, "Icohp_bndg_orb_sum_avg"] = 0 if len(icohp_sum_orb_bndg) == 0 else np.mean(icohp_sum_orb_bndg)
            df.loc[ids, "Icohp_bndg_orb_sum_max"] = 0 if len(icohp_sum_orb_bndg) == 0 else np.max(icohp_sum_orb_bndg)
            df.loc[ids, "Icohp_bndg_orb_sum_min"] = 0 if len(icohp_sum_orb_bndg) == 0 else np.min(icohp_sum_orb_bndg)
            df.loc[ids, "Icohp_bndg_orb_sum_std"] = 0 if len(icohp_sum_orb_bndg) == 0 else np.std(icohp_sum_orb_bndg)

            df.loc[ids, "bonding_orb_perc_avg"] = 0 if len(bond_orb) == 0 else np.mean(bond_orb)
            df.loc[ids, "bonding_orb_perc_max"] = 0 if len(bond_orb) == 0 else np.max(bond_orb)
            df.loc[ids, "bonding_orb_perc_min"] = 0 if len(bond_orb) == 0 else np.min(bond_orb)
            df.loc[ids, "bonding_orb_perc_std"] = 0 if len(bond_orb) == 0 else np.std(bond_orb)

            # anti-bonding orbital
            df.loc[ids, "Icohp_antibndg_orb_mean_avg"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.mean(icohp_mean_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_mean_max"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.max(icohp_mean_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_mean_min"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.min(icohp_mean_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_mean_std"] = (
                0 if len(icohp_mean_orb_antibndg) == 0 else np.std(icohp_mean_orb_antibndg)
            )

            df.loc[ids, "Icohp_antibndg_orb_sum_avg"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.mean(icohp_sum_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_sum_max"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.max(icohp_sum_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_sum_min"] = (
                0 if len(icohp_sum_orb_antibndg) == 0 else np.min(icohp_sum_orb_antibndg)
            )
            df.loc[ids, "Icohp_antibndg_orb_sum_std"] = (
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
    def get_lobsterpy_cba_dict(path_to_lobster_calc: str | Path, bonds: str, orbital_resolved: bool) -> dict:
        """
        Generate a Python dictionary object using the Analysis class with condensed bonding analysis data.

        :param path_to_lobster_calc: path to lobsterpy lightweight json file
        :param bonds: "all" or "cation-anion" bonds
        :param orbital_resolved:  bool indicating whether analysis is performed orbital wise

        Returns:
            Returns a dictionary with lobster summarized bonding analysis data

        """
        file_paths = get_file_paths(
            path_to_lobster_calc=path_to_lobster_calc, requested_files=["structure", "cohpcar", "icohplist", "charge"]
        )

        which_bonds = bonds.replace("-", "_")
        bond_type = f"{which_bonds}_bonds"

        try:
            analyse = Analysis(
                path_to_poscar=str(file_paths.get("structure")),
                path_to_icohplist=str(file_paths.get("icohplist")),
                path_to_cohpcar=str(file_paths.get("cohpcar")),
                path_to_charge=str(file_paths.get("charge")),
                summed_spins=False,  # we will always use spin polarization here
                cutoff_icohp=0.10,
                which_bonds=which_bonds,
                orbital_resolved=orbital_resolved,
            )

            data = {bond_type: {"lobsterpy_data": analyse.condensed_bonding_analysis}}
        except ValueError:
            data = {bond_type: {"lobsterpy_data": {}}}

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
                " Will set Madelung Energies for crystal structure values to NaN"
            )
            madelung_energies = {
                "Mulliken": np.nan,
                "Loewdin": np.nan,
                "Ewald_splitting": np.nan,
            }

            data["madelung_energies"] = madelung_energies

        return data


class CoxxFingerprint(NamedTuple):
    """
    Represents a Coxx fingerprint.

    This named tuple is used to store information related to a Coxx fingerprint, which
    includes energies, Coxx values, fingerprint type, spin type, number of bins, and bin width.

    :param energies: The energy values associated with the Coxx fingerprint.
    :param coxx: The Coxx values corresponding to each energy.
    :param fp_type: The type of the Coxx fingerprint.
    :param spin_type: The spin type associated with the fingerprint.
    :param n_bins: The number of bins used in the Coxx fingerprint.
    :param bin_width: The width of each bin in the Coxx fingerprint.
    """

    energies: np.ndarray
    coxx: np.ndarray
    fp_type: str
    spin_type: str
    n_bins: int
    bin_width: float


class FeaturizeCOXX:
    """
    Class to featurize COHPCAR, COBICAR or COOPCAR data.

    :param path_to_coxxcar: path to COXXCAR.lobster (e.g., `COXXCAR.lobster`)
    :param path_to_icoxxlist: path to ICOXXLIST.lobster (e.g., `ICOXXLIST.lobster`)
    :param path_to_structure: path to structure file (e.g., `POSCAR`)
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
        feature_type: str,
        e_range: list[float] = [-10.0, 0.0],
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
        """
        Featurize COHPCAR, COBICAR or COOPCAR data.

        :param path_to_coxxcar: path to COXXCAR.lobster (e.g., `COXXCAR.lobster`)
        :param path_to_icoxxlist: path to ICOXXLIST.lobster (e.g., `ICOXXLIST.lobster`)
        :param path_to_structure: path to structure file (e.g., `POSCAR`)
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
        spin_type: str = "summed",
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
        feature_type: str,
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

        return np.trapz(p**n * coxx, x=energies) / np.trapz(coxx, x=energies)

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


class FeaturizeCharges:
    """
    Class to compute ionicity from CHARGE.lobster data.

    :param path_to_structure: path to POSCAR
    :param path_to_charge: path to CHARGE.lobster (e.g., `CHARGE.lobster`)
    :param charge_type: set charge type used for computing ionicity.
        Possible options are `Mulliken` or `Loewdin`
    """

    def __init__(
        self,
        path_to_structure: str | Path,
        path_to_charge: str | Path,
        charge_type: str,
    ):
        """
        Compute the Ionicity of the structure from CHARGE.lobster data.

        :param path_to_structure: path to POSCAR
        :param path_to_charge: path to CHARGE.lobster (e.g., `CHARGE.lobster`)
        :param charge_type: set charge type used for computing ionicity.
            Possible options are `Mulliken` or `Loewdin`
        """
        self.path_to_structure = path_to_structure
        self.path_to_charge = path_to_charge
        self.charge_type = charge_type

    def _calc_ionicity(self) -> float:
        r"""
        Calculate ionicity of the crystal structure based on quantum chemical charges.

        \\[I_{\text{{Charges}}} = \frac{1}{{N_{\text{{Atoms}}}}}\\sum_{i=1}^{N_{\text{{Atoms}}} }
        \\left(\frac{q_i}{v_{\text{{eff}},i}}\right)\\]

        Returns:
            Ionicity of the structure

        """
        chargeobj = Charge(filename=self.path_to_charge)
        structure = Structure.from_file(self.path_to_structure)

        if self.charge_type.lower() not in ["mulliken", "loewdin"]:
            raise ValueError("Please check the requested charge_type. Possible options are `mulliken` or `loewdin`")

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

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        Return a pandas dataframe with computed ionicity as columns.

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.

        Returns:
            Returns a pandas dataframe with ionicity

        """
        if ids:
            df = pd.DataFrame(index=[ids])
        else:
            ids = Path(self.path_to_charge).parent.name
            df = pd.DataFrame(index=[ids])

        if self.charge_type.lower() == "mulliken":
            df.loc[ids, "Ionicity_Mull"] = self._calc_ionicity()
        else:
            df.loc[ids, "Ionicity_Loew"] = self._calc_ionicity()

        return df


class DosFingerprint(NamedTuple):
    """
    Represents a Density of States (DOS) fingerprint.

    This named tuple is used to store information related to the Density of States (DOS)
    in a material. It includes the energies, densities, type, number of bins, and bin width.

    :param energies: The energy values associated with the DOS.
    :param densities: The corresponding density values for each energy.
    :param type: The type of DOS fingerprint.
    :param n_bins: The number of bins used in the fingerprint.
    :param bin_width: The width of each bin in the DOS fingerprint.
    """

    energies: np.ndarray
    densities: np.ndarray
    type: str
    n_bins: int
    bin_width: float


class FeaturizeDoscar:
    """
    Class to compute DOS moments and fingerprints from DOSCAR.lobster / DOSCAR.LSO.lobster.

    :param path_to_structure: path to POSCAR
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

        :param path_to_structure: path to POSCAR
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
        fp_type: str = "summed_pdos",
        binning: bool = True,
        n_bins: int = 256,
        normalize: bool = True,
    ) -> pd.DataFrame:
        """
        Generate a dataframe consisting of DOS fingerprint (fp).

        :param ids: set index name in the pandas dataframe. Default is None.
            When None, LOBSTER calc directory name is used as index name.
        :param fp_type: Specify fingerprint type to compute, can accept `s/p/d/f/summed_pdos`
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

        fp = self.dos.get_dos_fp(
            type=fp_type,
            normalize=normalize,
            n_bins=n_bins,
            binning=binning,
            max_e=self.e_range[-1] if self.e_range is not None else None,
            min_e=self.e_range[0] if self.e_range is not None else None,
        )._asdict()

        df.loc[ids, "DOS_FP"] = DosFingerprint(**fp)

        return df
