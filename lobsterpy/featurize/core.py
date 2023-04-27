# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
This module defines functions to featurize Lobster lighweight jsons
"""

import gzip
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Mapping, NamedTuple, List
from collections import namedtuple
from pymatgen.core.structure import Structure
from mendeleev import element
from pymatgen.io.lobster import Charge, Lobsterout, Icohplist
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from scipy.integrate import trapezoid


class FeaturizeLobsterpy:
    """
    class to featurize lobsterpy data

    Args:
        path_to_json: path to lobster lightweight json
        bonds: "all_bonds" or "cation_anion_bonds"
    Attributes:
        get_df: returns a pandas dataframe with relevant icohp statistical data as columns from
        lobsterpy automatic bonding analysis

    """

    def __init__(
        self,
        path_to_json: str,
        bonds: str = "all_bonds",
    ):
        self.path_to_json = path_to_json
        self.bonds = bonds

    def get_df(self, ids: str) -> pd.DataFrame:
        """
        This function featurizes LobsterPy condensed bonding analysis data from
        lobster lightweight json.gz files

        Returns:
            Returns a pandas dataframe with lobsterpy icohp statistics

        """
        # define a pandas dataframe
        df = pd.DataFrame(index=[ids])

        # read the lightweight lobster json files using read_lobster_lightweight_json function
        data = self._read_lobster_lightweight_json()

        icohp_mean = []
        icohp_sum = []
        bond = []
        antibond = []
        # extract lobsterpy icohp related data for bond type specified
        if data[self.bonds]:
            for k, v in data[self.bonds]["lobsterpy_data"]["sites"].items():
                for k1, v1 in v["bonds"].items():
                    icohp_mean.append(float(v1["ICOHP_mean"]))
                    icohp_sum.append(float(v1["ICOHP_sum"]))
                    bond.append(v1["bonding"]["perc"])
                    antibond.append(v1["antibonding"]["perc"])
        else:
            raise Exception(
                "There exist no data for {} for this structure. " 'Please switch to "all_bonds" mode'.format(self.bonds)
            )

        # add stats data as columns to the dataframe
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

    def _read_lobster_lightweight_json(self):
        """
        This function reads loads the lightweight json.gz files and returns a python dictionary object
        with lobster summmarized bonding analysis data.

        Returns:
            Returns a dictionary with lobster summmarized bonding analysis data

        """
        with gzip.open(self.path_to_json, "rb") as f:
            data = json.loads(f.read().decode("utf-8"))

        lobster_data = {}
        for item in data:
            lobster_data.update(item)

        return lobster_data


class FeaturizeCOXX:

    """
    class to featurize  COHPCAR/COBICAR/COOPCAR data as fingerprint

    Args:
        path_to_coxxcar: path to COXXCAR.lobster (e.g., "COXXCAR.lobster")
        path_to_icoxxlist : path to ICOXXLIST.lobster (e.g., "ICOXXLIST.lobster")
        path_to_structure : path to structure file (e.g., "POSCAR")
        feature_type: set the feature type for moment features and fingerprints. Possible options are "bonding", "antibonding" or "overall"
        are_cobis : bool indicating if file contains COBI/ICOBI data
        are_coops : bool indicating if file contains COOP/ICOOP data
        include_wicoxx : if True will include weighted ICOXX, effective co-ordination number
        in the featurized dataframe
        include_moment_features : if True will include ICOXX center,width, skew and kurtosis features in the  featurized dataframe
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
        are_cobis: bool = False,
        are_coops: bool = False,
        include_wicoxx: bool = False,
        include_moment_features: bool = False,
        e_range: list = [-5, 0],
    ):
        self.path_to_coxxcar = path_to_coxxcar
        self.path_to_icoxxlist = path_to_icoxxlist
        self.path_to_structure = path_to_structure
        self.feature_type = feature_type
        self.include_wicoxx = include_wicoxx
        self.include_moment_features = include_moment_features
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
        )


    def get_coxx_fingerprint_df(
        self,
        ids: str,
        fp_type: str,
        label_list: list,
        coxxcar_obj: CompleteCohp,
        spin_type: str = "summed",
        binning: bool = True,
        min_e: float = -5,
        max_e: float = 0,
        n_bins: int = 50,
        normalize: bool = True,
    ):
        """
        Generates the COHP fingerprint

        Args:
            ids: sets index of pandas dataframe
            spin_type: Specify spin type. Can accept '{summed/up/down}'
            (default is summed)
            binning: If true cohps will be binned
            min_e: The minimum mode energy to include in the fingerprint (default is None)
            max_e: The maximum mode energy to include in the fingerprint (default is None)
            n_bins: Number of bins to be used in the fingerprint (default is 256)
            normalize: If true, normalizes the area under fp to equal to 1 (default is True)
            fp_type: Specify fingerprint type needed can accept '{bonding/antibonding/overall}'
            (default is antibonding)
            label_list: Specify bond lables as a list for which cohp fingerprints are needed
            coxxcar_obj: coxxcar object

        Raises:
            Exception: If spin_type is not one of the accepted values {summed/up/down}.
            ValueError: If fp_type is not one of the accepted values {bonding/antibonding/overall}.

        Returns:
            Fingerprint(namedtuple) : The COHP fingerprint
            of format (energies, coxx, fp_type, spin_type, n_bins, bin_width)
        """

        coxx_fingerprint = namedtuple("coxx_fingerprint", "energies coxx fp_type spin_type n_bins bin_width")
        energies = coxxcar_obj.energies - coxxcar_obj.efermi  # here substraction of efermi has impact on fp.

        if max_e is None:
            max_e = np.max(energies)

        if min_e is None:
            min_e = np.min(energies)

        if label_list:
            coxxcar_obj = coxxcar_obj.get_summed_cohp_by_label_list(label_list).get_cohp()
        else:
            coxxcar_obj = coxxcar_obj.get_cohp()

        if spin_type == "up":
            coxx_all = coxxcar_obj[Spin.up]
        elif spin_type == "down":
            coxx_all = coxxcar_obj[Spin.down]
        elif spin_type == "summed":
            coxx_all = coxxcar_obj[Spin.up] + coxxcar_obj[Spin.down]
        else:
            raise Exception(
                """Check the spin_type argument,
                                 Possible options are summed/up/down"""
            )
        coxx_dict = dict()

        coxx_dict["bonding"] = np.array([scohp if scohp <= 0 else 0 for scohp in coxx_all])
        coxx_dict["antibonding"] = np.array([scohp if scohp >= 0 else 0 for scohp in coxx_all])
        coxx_dict["overall"] = coxx_all

        try:
            coxxs = coxx_dict[fp_type]
            if len(energies) < n_bins:
                inds = np.where((energies >= min_e) & (energies <= max_e))
                return coxx_fingerprint(
                    energies[inds],
                    coxxs[inds],
                    fp_type,
                    spin_type,
                    len(energies),
                    np.diff(energies)[0],
                )

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
                fp_type,
                spin_type,
                n_bins,
                bin_width,
            )

            df_temp = pd.DataFrame(index=[ids], columns=["COHP_FP"])
            df_temp.at[ids, "COHP_FP"] = fp.coxx

            df = pd.DataFrame(df_temp["COHP_FP"].tolist())

            return df

        except KeyError:
            raise ValueError(
                """Please recheck fp_type requested argument.Possible options are bonding/antibonding/overall
                """
            )

    def _calculate_wicoxx_ein(self):
        """
        Method to calculate weighted icoxx (xx ==hp,op,bi) and ein of crystal structure based on work by Peter Müller

        {:[ bar(COHP)(E)=sum _(i)(w_(i)*COHP_(i)(E))],[w_(i)=(ICOHP_(i))/(ICOHP_("total "))]:}

        w_(i)=(ICOHP_(i))/(ICOHP_("total "))

        bar(ICOHP)=int_(-oo)^(epsi_(F)) bar(COHP)(E)dE

        Returns:
            Percent bonding, Percent anti-bonding, weighted icoxx, effective interaction number

        """
        icoxx_dict = self.icoxxlist.icohpcollection.as_dict()
        list_labels = icoxx_dict["list_labels"]
        # Compute sum of icohps
        icoxx_total = self.icoxxlist.icohpcollection.get_summed_icohp_by_label_list(list_labels)

        summed_weighted_coxx = []

        for lab in list_labels:
            for k, v in self.completecoxx.get_cohp_by_label("{}".format(lab), summed_spin_channels=True).cohp.items():
                coxx = v
            icoxx = self.icoxxlist.icohpcollection.get_icohp_by_label(lab, summed_spin_channels=True)
            weight = icoxx / icoxx_total  # calculate the weights based on icohp contri to total icohp of the structure
            weighted_coxx = weight * coxx
            summed_weighted_coxx.append(weighted_coxx)

        summed = np.sum(summed_weighted_coxx, axis=0)
        # below fermi
        indices = self.completecoxx.energies <= self.completecoxx.efermi
        en_bf = self.completecoxx.energies[indices]
        coxx_bf = summed[indices]
        # percent bonding-anitbonding
        if not self.icoxxlist.are_cobis and not self.icoxxlist.are_coops:
            bonding_indices = coxx_bf <= 0
            antibonding_indices = coxx_bf >= 0
            bnd = abs(trapezoid(en_bf[bonding_indices], coxx_bf[bonding_indices]))
            antibnd = abs(trapezoid(en_bf[antibonding_indices], coxx_bf[antibonding_indices]))
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100
        elif self.icoxxlist.are_cobis or self.icoxxlist.are_coops:
            bonding_indices = coxx_bf >= 0
            antibonding_indices = coxx_bf <= 0
            bnd = abs(trapezoid(coxx_bf[bonding_indices], en_bf[bonding_indices]))
            antibnd = abs(trapezoid(coxx_bf[antibonding_indices], en_bf[antibonding_indices]))
            per_bnd = (bnd / (bnd + antibnd)) * 100
            per_antibnd = (antibnd / (bnd + antibnd)) * 100

        w_icoxx = trapezoid(coxx_bf, en_bf)

        ein = (icoxx_total / w_icoxx) * (2 / self.completecoxx.structure.num_sites)  # calc effective interaction number

        return per_bnd, per_antibnd, w_icoxx, ein

    def _calc_moment_features(self):
        coxxcar = self.completecoxx.get_cohp()
        coxx_all = coxxcar[Spin.up] + coxxcar[Spin.down]
        energies = self.completecoxx.energies - self.completecoxx.efermi

        coxx_dict = {}
        if not self.are_cobis and not self.are_coops:
            coxx_dict["bonding"] = np.array([scohp if scohp <= 0 else 0 for scohp in coxx_all])
            coxx_dict["antibonding"] = np.array([scohp if scohp >= 0 else 0 for scohp in coxx_all])
            coxx_dict["overall"] = coxx_all
        else:
            coxx_dict["antibonding"] = np.array([scohp if scohp <= 0 else 0 for scohp in coxx_all])
            coxx_dict["bonding"] = np.array([scohp if scohp >= 0 else 0 for scohp in coxx_all])
            coxx_dict["overall"] = coxx_all

        band_center = self._get_band_center(
            cohp=coxx_dict[self.feature_type],
            energies=energies,
            e_range=self.e_range,
        )
        band_width = self._get_n_moment(
            n=2,
            cohp=coxx_dict[self.feature_type],
            energies=energies,
            e_range=self.e_range,
        )
        band_skew = self._get_n_moment(
            n=3,
            cohp=coxx_dict[self.feature_type],
            energies=energies,
            e_range=self.e_range,
        )
        band_kurt = self._get_n_moment(
            n=4,
            cohp=coxx_dict[self.feature_type],
            energies=energies,
            e_range=self.e_range,
        )

        return band_center, band_width, band_skew, band_kurt

    def _get_band_center(self, cohp: List[float], energies: List[float], e_range: List[int]) -> float:
        """
        Get the band width, defined as the first moment of the orbital resolved COHP
        Args:
            cohp: Orbital COHP list
            energies: energies corresponding orbital orbital interactions COHP
        Returns:
            Orbital-Orbital interaction band center in eV
        """
        band_center = self._get_n_moment(n=1, cohp=cohp, energies=energies, e_range=e_range, center=False)

        return band_center

    def _get_n_moment(self, n: int, cohp: list, energies: list, e_range: list, center: bool = True) -> float:
        """

        Get the nth moment of COXX

        Args:
            n: The order for the moment
            cohp: COXX list
            energies: energies
            center: Take moments with respect to the COXX center

        Returns:
            COXX nth moment in eV
        """
        if e_range:
            cohp = cohp[(energies >= self.e_range[0]) & (energies <= self.e_range[1])]
            energies = energies[(energies >= self.e_range[0]) & (energies <= self.e_range[1])]

        if center:
            band_center = self._get_band_center(cohp=cohp, energies=energies, e_range=self.e_range)
            p = energies - band_center
        else:
            p = energies

        nth_moment = np.trapz(p**n * cohp, x=energies) / np.trapz(cohp, x=energies)

        return nth_moment

    def get_df(self, ids: str):
        """
        This function returns a pandas dataframe with weighted ICOXX, effective interaction number
        and moment features (center, width, skewness and kurtosis) of COXX in selected energy range

        Returns:
            Returns a pandas dataframe with cohp/cobi/coop related features

        """

        df = pd.DataFrame(index=[ids])

        if self.include_wicoxx:
            (
                per_bnd_xx,
                per_antibnd_xx,
                w_icoxx,
                ein_xx,
            ) = self._calculate_wicoxx_ein()

            if self.are_coops:
                bc, bw, bs, bk = self._calc_moment_features()
                df.loc[ids, "bnd_ICOOP"] = per_bnd_xx
                df.loc[ids, "antibnd_ICOOP"] = per_antibnd_xx
                df.loc[ids, "w_ICOOP"] = w_icoxx
                df.loc[ids, "EIN_ICOOP"] = ein_xx
                df.loc[ids, "band_center_COOP"] = bc
                df.loc[ids, "band_width_COOP"] = bw
                df.loc[ids, "band_skewness_COOP"] = bs
                df.loc[ids, "band_kurtosis_COOP"] = bk
            elif self.are_cobis:
                bc, bw, bs, bk = self._calc_moment_features()
                df.loc[ids, "bnd_ICOOP"] = per_bnd_xx
                df.loc[ids, "bnd_ICOBI"] = per_bnd_xx
                df.loc[ids, "antibnd_ICOBI"] = per_antibnd_xx
                df.loc[ids, "w_ICOBI"] = w_icoxx
                df.loc[ids, "EIN_ICOBI"] = ein_xx
                df.loc[ids, "band_center_COBi"] = bc
                df.loc[ids, "band_width_COBI"] = bw
                df.loc[ids, "band_skewness_COBI"] = bs
                df.loc[ids, "band_kurtosis_COBI"] = bk
            else:
                bc, bw, bs, bk = self._calc_moment_features()
                df.loc[ids, "bnd_ICOHP"] = per_bnd_xx
                df.loc[ids, "antibnd_ICOHP"] = per_antibnd_xx
                df.loc[ids, "w_ICOHP"] = w_icoxx
                df.loc[ids, "EIN_ICOHP"] = ein_xx
                df.loc[ids, "band_center_COHP"] = bc
                df.loc[ids, "band_width_COHP"] = bw
                df.loc[ids, "band_skewness_COHP"] = bs
                df.loc[ids, "band_kurtosis_COHP"] = bk

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
        path_to_strucutre: str,
        path_to_charge: str,
        charge_type: str,
    ):
        self.path_to_structure = path_to_strucutre
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

        if self.charge_type.lower() == "mulliken":
            charges = chargeobj.Mulliken
        elif self.charge_type.lower() == "loewdin":
            charges = chargeobj.Loewdin
        else:
            raise ValueError("Please check the requested charge_type, " 'Possible options are "Mulliken" or "Loewdin"')

        ch_veff = []
        for i, j in enumerate(charges):
            if (
                j > 0
                and not structure.species[i].is_transition_metal
                and (not structure.species[i].is_actinoid and not structure.species[i].is_lanthanoid)
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 0)
                ch_veff.append(val)

            elif (
                j > 0
                and not structure.species[i].is_transition_metal
                and (not structure.species[i].is_actinoid and not structure.species[i].is_lanthanoid)
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 8)
                ch_veff.append(val)

            elif (
                j < 0
                and not structure.species[i].is_transition_metal
                and (not structure.species[i].is_actinoid and not structure.species[i].is_lanthanoid)
            ):
                valence_elec = element(structure.species[i].value)
                val = j / (valence_elec.nvalence() - 8)
                ch_veff.append(val)

            elif j > 0 and structure.species[i].is_transition_metal:
                val = j / (structure.species[i].max_oxidation_state - 0)
                ch_veff.append(val)

            elif j < 0 and structure.species[i].is_transition_metal:
                val = j / (structure.species[i].max_oxidation_state - 8)
                ch_veff.append(val)

            elif j > 0 and structure.species[i].is_lanthanoid:
                val = j / (structure.species[i].max_oxidation_state - 0)
                ch_veff.append(val)

            elif j < 0 and structure.species[i].is_actinoid:
                val = j / (structure.species[i].max_oxidation_state - 8)
                ch_veff.append(val)

        ionicity = sum(ch_veff) / structure.num_sites

        return ionicity

    def get_df(self, ids):
        """
        This function returns a pandas dataframe with computed ionicity as column

        Returns:
            Returns a pandas dataframe with ionicity

        """

        df = pd.DataFrame(index=[ids])

        if self.charge_type == "mulliken":
            df.loc[ids, "Ionicity_Mull"] = self._calc_ionicity()
        else:
            df.loc[ids, "Ionicity_Loew"] = self._calc_ionicity()

        return df
