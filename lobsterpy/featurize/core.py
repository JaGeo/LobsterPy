# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines classes to featurize Lobster data ready to be used for ML studies
"""

from __future__ import annotations
import gzip
import json
import os
import warnings
from pathlib import Path
from typing import List, Tuple
from collections import namedtuple
import numpy as np
import numpy.typing as npt
import pandas as pd
from mendeleev import element
from pymatgen.core.structure import Structure
from pymatgen.io.lobster import Charge, Icohplist, MadelungEnergies
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from scipy.integrate import trapezoid
from lobsterpy.cohp.analyze import Analysis

warnings.filterwarnings("ignore")


class FeaturizeLobsterpy:
    """
    class to featurize lobsterpy data

    Args:
        path_to_lobster_calc: path containing lobster calc outputs
        path_to_json: path to lobster lightweight json
        bonds: "all" or "cation-anion" bonds
    Attributes:
        get_df: returns a pandas dataframe with relevant icohp statistical data as columns from
        lobsterpy automatic bonding analysis

    """

    def __init__(
        self,
        path_to_lobster_calc: str | None = None,
        path_to_json: str | None = None,
        bonds: str = "all",
    ):
        self.path_to_json = path_to_json
        self.path_to_lobster_calc = path_to_lobster_calc
        self.bonds = bonds

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        This function featurizes LobsterPy condensed bonding analysis data from
        lobster lightweight json.gz files

        Returns:
            Returns a pandas dataframe with lobsterpy icohp statistics

        """
        if self.path_to_json and not self.path_to_lobster_calc:
            # read the lightweight lobster json files using read_lobster_lightweight_json method
            data = FeaturizeLobsterpy.read_lobster_lightweight_json(
                path_to_json=self.path_to_json
            )
            if not ids:
                ids = Path(self.path_to_json).name.split(".")[0]

        elif self.path_to_lobster_calc and not self.path_to_json:
            # get lobsterpy condensed bonding analysis data using get_lobsterpy_cba_dict method
            data = FeaturizeLobsterpy.get_lobsterpy_cba_dict(
                path_to_lobster_calc=self.path_to_lobster_calc, bonds=self.bonds
            )

            if not ids:
                ids = Path(self.path_to_lobster_calc).name

        else:
            raise ValueError(
                "Please provide either path to lightweight lobster jsons or path to lobster calc"
            )
        # define a pandas dataframe
        df = pd.DataFrame(index=[ids])

        icohp_mean = []
        icohp_sum = []
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

        if (
            not data[bond_type]["lobsterpy_data"]
            or not data[bond_type]["lobsterpy_data"]["sites"]
        ):
            raise Exception(
                "No {} bonds detected for {} structure. "
                "Please switch to ´all´ bonds mode".format(self.bonds, ids)
            )

        for k, v in data[bond_type]["lobsterpy_data"]["sites"].items():
            if v["bonds"]:
                for k1, v1 in v["bonds"].items():
                    icohp_mean.append(float(v1["ICOHP_mean"]))
                    icohp_sum.append(float(v1["ICOHP_sum"]))
                    bond.append(v1["bonding"]["perc"])
                    antibond.append(v1["antibonding"]["perc"])

        # add ICOHP stats data (mean, min, max, standard deviation) as columns to the dataframe
        df.loc[ids, "Icohp_mean_avg"] = np.mean(icohp_mean)
        df.loc[ids, "Icohp_mean_max"] = np.max(icohp_mean)
        df.loc[ids, "Icohp_mean_min"] = np.min(icohp_mean)
        df.loc[ids, "Icohp_mean_std"] = np.std(icohp_mean)

        df.loc[ids, "Icohp_sum_avg"] = np.mean(icohp_sum)
        df.loc[ids, "Icohp_sum_max"] = np.max(icohp_sum)
        df.loc[ids, "Icohp_sum_min"] = np.min(icohp_sum)
        df.loc[ids, "Icohp_sum_std"] = np.std(icohp_sum)

        df.loc[ids, "bonding_perc_avg"] = np.mean(bond)
        df.loc[ids, "bonding_perc_max"] = np.max(bond)
        df.loc[ids, "bonding_perc_min"] = np.min(bond)
        df.loc[ids, "bonding_perc_std"] = np.std(bond)

        df.loc[ids, "antibonding_perc_avg"] = np.mean(antibond)
        df.loc[ids, "antibonding_perc_min"] = np.min(antibond)
        df.loc[ids, "antibonding_perc_max"] = np.max(antibond)
        df.loc[ids, "antibonding_perc_std"] = np.std(antibond)

        # add madelung energies for the structure
        df.loc[ids, "Madelung_Mull"] = data["madelung_energies"]["Mulliken"]
        df.loc[ids, "Madelung_Loew"] = data["madelung_energies"]["Loewdin"]

        return df

    @staticmethod
    def read_lobster_lightweight_json(path_to_json: str) -> dict:
        """
        This method reads loads the lightweight json.gz files and returns a python dictionary object
        with lobster summmarized bonding analysis data.

        Args:
            path_to_json: path to lobsterpy lightweight json file

        Returns:
            Returns a dictionary with lobster summmarized bonding analysis data

        """
        with gzip.open(str(path_to_json), "rb") as f:
            data = json.loads(f.read().decode("utf-8"))

        lobster_data = {}
        for item in data:
            lobster_data.update(item)

        return lobster_data

    @staticmethod
    def get_lobsterpy_cba_dict(path_to_lobster_calc: str, bonds: str) -> dict:
        """
        This function uses lobsterpy.cohp.analyze.Analysis class to generate a python dictionary object
        with lobster summmarized bonding analysis data.

        Args:
            path_to_lobster_calc: path to lobsterpy lightweight json file
            bonds: "all" or "cation-anion" bonds

        Returns:
            Returns a dictionary with lobster summmarized bonding analysis data

        """
        dir_name = Path(str(path_to_lobster_calc))

        # check if files are compressed (.gz) and update file paths
        req_files_lobsterpy = {
            "structure_path": "POSCAR",
            "cohpcar_path": "COHPCAR.lobster",
            "icohplist_path": "ICOHPLIST.lobster",
            "charge_path": "CHARGE.lobster",
        }

        for file, default_value in req_files_lobsterpy.items():
            file_path = dir_name / default_value
            req_files_lobsterpy[file] = file_path  # type: ignore
            if not file_path.exists():
                gz_file_path = file_path.with_name(file_path.name + ".gz")
                if gz_file_path.exists():
                    req_files_lobsterpy[file] = gz_file_path  # type: ignore
                else:
                    raise Exception(
                        "Path provided for Lobster calc directory seems incorrect."
                        "It does not contain COHPCAR.lobster, ICOHPLIST.lobster, POSCAR and "
                        "CHARGE.lobster files needed for automatic analysis using LobsterPy"
                    )

        cohpcar_path = req_files_lobsterpy.get("cohpcar_path")
        charge_path = req_files_lobsterpy.get("charge_path")
        structure_path = req_files_lobsterpy.get("structure_path")
        icohplist_path = req_files_lobsterpy.get("icohplist_path")

        which_bonds = bonds

        if which_bonds == "all":
            bond_type = "all_bonds"
        elif which_bonds == "cation-anion":
            bond_type = "cation_anion_bonds"

        try:
            analyse = Analysis(
                path_to_poscar=str(structure_path),
                path_to_icohplist=str(icohplist_path),
                path_to_cohpcar=str(cohpcar_path),
                path_to_charge=str(charge_path),
                summed_spins=False,  # we will always use spin polarization here
                cutoff_icohp=0.10,
                whichbonds=which_bonds,
            )

            data = {bond_type: {"lobsterpy_data": analyse.condensed_bonding_analysis}}
        except ValueError:
            data = {bond_type: {"lobsterpy_data": {}}}

        madelung_energies_path = dir_name / "MadelungEnergies.lobster"
        # check if .gz file exists and update Madelung Energies path
        if not madelung_energies_path.exists():
            gz_file_path = madelung_energies_path.with_name(
                madelung_energies_path.name + ".gz"
            )
            if gz_file_path.exists():
                madelung_energies_path = gz_file_path

        if madelung_energies_path.exists():
            madelung_obj = MadelungEnergies(filename=str(madelung_energies_path))

            madelung_energies = {
                "Mulliken": madelung_obj.madelungenergies_Mulliken,
                "Loewdin": madelung_obj.madelungenergies_Loewdin,
                "Ewald_splitting": madelung_obj.ewald_splitting,
            }
            data["madelung_energies"] = madelung_energies

        else:
            warnings.warn(
                "MadelungEnergies.lobster file not found in Lobster calc directory provided"
                " Will set Madelung Engeries for crystal structure values to NaN"
            )
            madelung_energies = {
                "Mulliken": np.nan,
                "Loewdin": np.nan,
                "Ewald_splitting": np.nan,
            }

            data["madelung_energies"] = madelung_energies

        return data


coxx_fingerprint = namedtuple(
    "coxx_fingerprint", "energies coxx fp_type spin_type n_bins bin_width"
)


class FeaturizeCOXX:
    """
    class to generate features from COHPCAR/COBICAR/COOPCAR data

    Args:
        path_to_coxxcar: path to COXXCAR.lobster (e.g., "COXXCAR.lobster")
        path_to_icoxxlist : path to ICOXXLIST.lobster (e.g., "ICOXXLIST.lobster")
        path_to_structure : path to structure file (e.g., "POSCAR")
        feature_type: set the feature type for moment features and fingerprints.
        Possible options are "bonding", "antibonding" or "overall"
        are_cobis : bool indicating if file contains COBI/ICOBI data
        are_coops : bool indicating if file contains COOP/ICOOP data
        e_range : range of energy relative to fermi for which moment features needs to be computed

    Attributes:
        get_df: pandas dataframe
        get_coxx_fingerprint_df: pandas dataframe

    """

    def __init__(
        self,
        path_to_coxxcar: str,
        path_to_icoxxlist: str,
        path_to_structure: str,
        feature_type: str,
        e_range: List[float] = [-10.0, 0.0],
        are_cobis: bool = False,
        are_coops: bool = False,
    ):
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
        label_list: List[str] | None = None,
        per_bond: bool = True,
        spin_type: str = "summed",
        binning: bool = True,
        n_bins: int = 56,
        normalize: bool = True,
    ) -> pd.DataFrame:
        """
        Generates the COXX fingerprint

        Args:
            ids: sets index of pandas dataframe
            spin_type: Specify spin type. Can accept '{summed/up/down}'
            (default is summed)
            binning: If true coxxs will be binned
            n_bins: Number of bins to be used in the fingerprint (default is 256)
            normalize: If true, normalizes the area under fp to equal to 1 (default is True)
            label_list: Specify bond lables as a list for which cohp fingerprints are needed
            per_bond: Will scale cohp values by number of bonds i.e length of label_list arg
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

        if label_list:
            divisor = len(label_list) if per_bond else 1
            coxxcar_obj = coxxcar_obj.get_summed_cohp_by_label_list(
                label_list, divisor=divisor
            ).get_cohp()
        else:
            coxxcar_obj = coxxcar_obj.get_cohp()

        if spin_type == "up":
            coxx_all = coxxcar_obj[Spin.up]
        elif spin_type == "down":
            if Spin.down in coxxcar_obj:
                coxx_all = coxxcar_obj[Spin.down]
            else:
                raise ValueError(
                    "LOBSTER calculation is non-spin polarized. Please switch spin_type to `up`"
                )
        elif spin_type == "summed":
            if Spin.down in coxxcar_obj:
                coxx_all = coxxcar_obj[Spin.up] + coxxcar_obj[Spin.down]
            else:
                coxx_all = coxxcar_obj[Spin.up]
        else:
            raise Exception(
                "Check the spin_type argument." "Possible options are summed/up/down"
            )

        coxx_dict = {}
        tol = 1e-6
        if not self.are_cobis and not self.are_coops:
            coxx_dict["bonding"] = np.array(
                [scohp if scohp <= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["antibonding"] = np.array(
                [scohp if scohp >= tol else 0 for scohp in coxx_all]
            )
        else:
            coxx_dict["antibonding"] = np.array(
                [scohp if scohp <= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["bonding"] = np.array(
                [scohp if scohp >= tol else 0 for scohp in coxx_all]
            )

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
                fp = coxx_fingerprint(
                    energies[inds],
                    coxxs[inds],
                    self.feature_type,
                    spin_type,
                    len(energies),
                    np.diff(energies)[0],
                )

                df.at[ids, "COXX_FP"] = fp

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

            fp = coxx_fingerprint(
                np.array([ener]),
                coxx_rebin_sc,
                self.feature_type,
                spin_type,
                n_bins,
                bin_width,
            )

            df.at[ids, "COXX_FP"] = fp

            return df

        except KeyError:
            raise ValueError(
                "Please recheck fp_type requested argument.Possible options are bonding/antibonding/overall"
            )

    def _calculate_wicoxx_ein(self) -> Tuple[float, float, float, float]:
        """
        Method to calculate weighted icoxx (xx ==hp,op,bi) and ein of crystal structure based on work by Peter Müller

        {:[ bar(COHP)(E)=sum _(i)(w_(i)*COHP_(i)(E))],[w_(i)=(ICOHP_(i))/(ICOHP_("total "))]:}

        w_(i)=(ICOHP_(i))/(ICOHP_("total "))

        bar(ICOHP)=int_(-oo)^(epsi_(F)) bar(COHP)(E)dE

        Returns:
            Percent bonding, Percent anti-bonding, weighted icoxx, effective interaction number

        """
        list_labels = list(self.icoxxlist.icohplist.keys())
        # Compute sum of icohps
        icoxx_total = self.icoxxlist.icohpcollection.get_summed_icohp_by_label_list(
            list_labels
        )

        summed_weighted_coxx = []

        for lab in list_labels:
            for k, v in self.completecoxx.get_cohp_by_label(
                "{}".format(lab), summed_spin_channels=True
            ).cohp.items():
                coxx = v
            icoxx = self.icoxxlist.icohpcollection.get_icohp_by_label(
                lab, summed_spin_channels=True
            )
            weight = (
                icoxx / icoxx_total
            )  # calculate the weights based on icohp contri to total icohp of the structure
            weighted_coxx = weight * coxx
            summed_weighted_coxx.append(weighted_coxx)

        summed = np.sum(summed_weighted_coxx, axis=0)
        # below fermi
        indices = self.completecoxx.energies <= self.completecoxx.efermi
        en_bf = self.completecoxx.energies[indices]
        coxx_bf = summed[indices]

        w_icoxx = trapezoid(coxx_bf, en_bf)

        ein = (icoxx_total / w_icoxx) * (
            2 / self.completecoxx.structure.num_sites
        )  # calc effective interaction number

        # percent bonding-anitbonding
        tol = 1e-6
        if not self.icoxxlist.are_cobis and not self.icoxxlist.are_coops:
            bonding_indices = coxx_bf <= tol
            antibonding_indices = coxx_bf >= tol
            bnd = abs(trapezoid(en_bf[bonding_indices], coxx_bf[bonding_indices]))
            antibnd = abs(
                trapezoid(en_bf[antibonding_indices], coxx_bf[antibonding_indices])
            )
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100

        elif self.icoxxlist.are_cobis or self.icoxxlist.are_coops:
            bonding_indices = coxx_bf >= tol
            antibonding_indices = coxx_bf <= tol
            bnd = abs(trapezoid(coxx_bf[bonding_indices], en_bf[bonding_indices]))
            antibnd = abs(
                trapezoid(coxx_bf[antibonding_indices], en_bf[antibonding_indices])
            )
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100

        return per_bnd, per_antibnd, w_icoxx, ein

    def _calc_moment_features(
        self, label_list: List[str] | None = None, per_bond=True
    ) -> Tuple[float, float, float, float]:
        """
        Wrapper method to calculate band center,width, skewness, and kurtosis of the COXX
        Args:
            label_list: List of bond labels
            per_bond: Will scale cohp values by number of bonds i.e length of label_list arg
            (Only affects when label_list is not None)

        Returns:
            coxx center,width, skewness, and kurtosis in eV
        """
        if label_list:
            divisor = len(label_list) if per_bond else 1
            coxxcar = self.completecoxx.get_summed_cohp_by_label_list(
                label_list, divisor=divisor
            ).get_cohp()

        else:
            coxxcar = self.completecoxx.get_cohp()

        if Spin.down in coxxcar:
            coxx_all = coxxcar[Spin.up] + coxxcar[Spin.down]
        else:
            coxx_all = coxxcar[Spin.up]

        energies = self.completecoxx.energies - self.completecoxx.efermi

        coxx_dict = {}
        tol = 1e-6
        if not self.are_cobis and not self.are_coops:
            coxx_dict["bonding"] = np.array(
                [scohp if scohp <= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["antibonding"] = np.array(
                [scohp if scohp >= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["overall"] = coxx_all
        else:
            coxx_dict["antibonding"] = np.array(
                [scohp if scohp <= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["bonding"] = np.array(
                [scohp if scohp >= tol else 0 for scohp in coxx_all]
            )
            coxx_dict["overall"] = coxx_all

        try:
            coxx_center = self._get_coxx_center(
                coxx=coxx_dict[self.feature_type],
                energies=energies,
                e_range=self.e_range,
            )
            coxx_width = np.sqrt(
                self._get_n_moment(
                    n=2,
                    coxx=coxx_dict[self.feature_type],
                    energies=energies,
                    e_range=self.e_range,
                )
            )
            coxx_skew = self._get_n_moment(
                n=3,
                coxx=coxx_dict[self.feature_type],
                energies=energies,
                e_range=self.e_range,
            ) / self._get_n_moment(
                n=2,
                coxx=coxx_dict[self.feature_type],
                energies=energies,
                e_range=self.e_range,
            ) ** (
                3 / 2
            )

            coxx_kurt = (
                self._get_n_moment(
                    n=4,
                    coxx=coxx_dict[self.feature_type],
                    energies=energies,
                    e_range=self.e_range,
                )
                / self._get_n_moment(
                    n=2,
                    coxx=coxx_dict[self.feature_type],
                    energies=energies,
                    e_range=self.e_range,
                )
                ** 2
            )
        except KeyError:
            raise ValueError(
                "Please recheck fp_type requested argument.Possible options are bonding/antibonding/overall"
            )

        return coxx_center, coxx_width, coxx_skew, coxx_kurt

    def _get_coxx_center(
        self,
        coxx: npt.NDArray[np.floating],
        energies: npt.NDArray[np.floating],
        e_range: List[float],
    ) -> float:
        """
        Get the band width, defined as the first moment of the COXX
        Args:
            coxx: COXX array
            energies: energies corresponding  COXX
            e_range: range of energy to compute coxx center

        Returns:
            coxx center in eV
        """
        coxx_center = self._get_n_moment(
            n=1, coxx=coxx, energies=energies, e_range=e_range, center=False
        )

        return coxx_center

    def _get_n_moment(
        self,
        n: float,
        coxx: npt.NDArray[np.floating],
        energies: npt.NDArray[np.floating],
        e_range: list[float] | None,
        center: bool = True,
    ) -> float:
        """

        Get the nth moment of COXX

        Args:
            n: The order for the moment
            coxx: COXX array
            energies: energies array
            center: Take moments with respect to the COXX center

        Returns:
            COXX nth moment in eV
        """
        if e_range:
            min_e = self.e_range[0]
            max_e = self.e_range[1]
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
            coxx_center = self._get_coxx_center(
                coxx=coxx, energies=energies, e_range=[min_e, max_e]
            )
            p = energies - coxx_center
        else:
            p = energies

        nth_moment = np.trapz(p**n * coxx, x=energies) / np.trapz(coxx, x=energies)  # type: ignore

        return nth_moment

    def get_summarized_coxx_df(
        self,
        ids: str | None = None,
        label_list: List[str] | None = None,
        per_bond=True,
    ) -> pd.DataFrame:
        """
        This function returns a pandas dataframe with weighted ICOXX, effective interaction number
        and moment features (center, width, skewness and kurtosis) of COXX in selected energy range

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

        if self.are_coops:
            cc, cw, cs, ck = self._calc_moment_features(
                label_list=label_list, per_bond=per_bond
            )
            df.loc[ids, "bnd_wICOOP"] = per_bnd_xx
            df.loc[ids, "antibnd_wICOOP"] = per_antibnd_xx
            df.loc[ids, "w_ICOOP"] = w_icoxx
            df.loc[ids, "EIN_ICOOP"] = ein_xx
            df.loc[ids, "center_COOP"] = cc
            df.loc[ids, "width_COOP"] = cw
            df.loc[ids, "skewness_COOP"] = cs
            df.loc[ids, "kurtosis_COOP"] = ck
        elif self.are_cobis:
            cc, cw, cs, ck = self._calc_moment_features(
                label_list=label_list, per_bond=per_bond
            )
            df.loc[ids, "bnd_wICOBI"] = per_bnd_xx
            df.loc[ids, "antibnd_wICOBI"] = per_antibnd_xx
            df.loc[ids, "w_ICOBI"] = w_icoxx
            df.loc[ids, "EIN_ICOBI"] = ein_xx
            df.loc[ids, "center_COBI"] = cc
            df.loc[ids, "width_COBI"] = cw
            df.loc[ids, "skewness_COBI"] = cs
            df.loc[ids, "kurtosis_COBI"] = ck
        else:
            cc, cw, cs, ck = self._calc_moment_features(
                label_list=label_list, per_bond=per_bond
            )
            df.loc[ids, "bnd_wICOHP"] = per_bnd_xx
            df.loc[ids, "antibnd_wICOHP"] = per_antibnd_xx
            df.loc[ids, "w_ICOHP"] = w_icoxx
            df.loc[ids, "EIN_ICOHP"] = ein_xx
            df.loc[ids, "center_COHP"] = cc
            df.loc[ids, "width_COHP"] = cw
            df.loc[ids, "skewness_COHP"] = cs
            df.loc[ids, "kurtosis_COHP"] = ck

        return df


class FeaturizeCharges:
    """
    class to compute ionicity from CHARGE.lobster data

    Args:
        path_to_strucutre: path to POSCAR
        path_to_charge : path to CHARGE.lobster (e.g., "CHARGE.lobster")
        charge_type : set charge type used for computing ionicity. Possible options are "Mulliken" or "Loewdin"

    Attributes:
        get_df: pandas dataframe

    """

    def __init__(
        self,
        path_to_structure: str,
        path_to_charge: str,
        charge_type: str,
    ):
        self.path_to_structure = path_to_structure
        self.path_to_charge = path_to_charge
        self.charge_type = charge_type

    def _calc_ionicity(self) -> float:
        """
        Method to calculate ionicity of crystal structure based on quantum chemical charges

        I_("Charges ")=(1)/(N_("Atoms "))sum _(i)^(N_("Atoms "))((q_(i))/(v_("eff ",i)))

        Returns:
            Ionicity of the structure

        """
        chargeobj = Charge(filename=self.path_to_charge)
        structure = Structure.from_file(self.path_to_structure)

        if self.charge_type.lower() not in ["mulliken", "loewdin"]:
            raise ValueError(
                "Please check the requested charge_type. "
                "Possible options are `Mulliken` or `Loewdin`"
            )

        ch_veff = []
        tol = 1e-6
        for i, j in enumerate(getattr(chargeobj, self.charge_type.capitalize())):
            if (
                j > tol
                and not structure.species[i].is_transition_metal
                and (
                    not structure.species[i].is_actinoid
                    and not structure.species[i].is_lanthanoid
                )
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 0)
                ch_veff.append(val)

            elif (
                j < tol
                and not structure.species[i].is_transition_metal
                and (
                    not structure.species[i].is_actinoid
                    and not structure.species[i].is_lanthanoid
                )
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 8)
                ch_veff.append(val)

            elif j > tol and (
                structure.species[i].is_transition_metal
                or structure.species[i].is_actinoid
                or structure.species[i].is_lanthanoid
            ):
                val = j / (structure.species[i].max_oxidation_state - 0)
                ch_veff.append(val)

            elif j < tol and (
                structure.species[i].is_transition_metal
                or structure.species[i].is_actinoid
                or structure.species[i].is_lanthanoid
            ):
                val = j / (structure.species[i].min_oxidation_state - 8)
                ch_veff.append(val)

        ionicity = sum(ch_veff) / structure.num_sites

        return ionicity

    def get_df(self, ids: str | None = None) -> pd.DataFrame:
        """
        This function returns a pandas dataframe with computed ionicity as column

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
