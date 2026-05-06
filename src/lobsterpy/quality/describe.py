"""This module contains functions to describe the quality of LOBSTER calculations."""


def get_calc_quality_description(quality_dict):
    """
    Generate a text description of the LOBSTER calculation quality.

    :param quality_dict: python dictionary from lobsterpy.analysis.get_lobster_calc_quality_summary
    """
    text_des = []

    for key, val in quality_dict.items():
        if key == "minimal_basis":
            if val:
                text_des.append("The LOBSTER calculation used minimal basis.")
            if not val:
                text_des.append(
                    "Consider rerunning the calculation with the minimum basis as well. Choosing a "
                    "larger basis set is only recommended if you see a significant improvement of "
                    "the charge spilling."
                )

        elif key == "charge_spilling":
            text_des.append(
                "The absolute and total charge spilling for the calculation is {} and {} %, respectively.".format(
                    quality_dict[key]["abs_charge_spilling"],
                    quality_dict[key]["abs_total_spilling"],
                )
            )
        elif key == "band_overlaps_analysis":
            if quality_dict[key]["file_exists"]:
                if quality_dict[key]["has_good_quality_maxDeviation"]:
                    text_des.append(
                        "The bandOverlaps.lobster file is generated during the LOBSTER run. This "
                        "indicates that the projected wave function is not completely orthonormalized; "
                        "however, the maximal deviation values observed compared to the identity matrix "
                        "is below the threshold of 0.1."
                    )
                else:
                    text_des.append(
                        "The bandOverlaps.lobster file is generated during the LOBSTER run. This "
                        "indicates that the projected wave function is not completely orthonormalized. "
                        "The maximal deviation value from the identity matrix is {}, and there are "
                        "{} percent k-points above the deviation threshold of 0.1. Please check the "
                        "results of other quality checks like dos comparisons, charges, "
                        "charge spillings before using the results for further "
                        "analysis.".format(
                            quality_dict[key]["max_deviation"],
                            quality_dict[key]["percent_kpoints_abv_limit"],
                        )
                    )
            else:
                text_des.append(
                    "The projected wave function is completely orthonormalized as no "
                    "bandOverlaps.lobster file is generated during the LOBSTER run."
                )

        elif key == "charge_comparisons":
            if val:
                for charge in ["mulliken", "loewdin"]:
                    if val[f"bva_{charge}_agree"]:
                        text_des.append(
                            f"The atomic charge signs from {charge.capitalize()} population analysis "
                            f"agree with the bond valence analysis."
                        )
                    if not val[f"bva_{charge}_agree"]:
                        text_des.append(
                            f"The atomic charge signs from {charge.capitalize()} population analysis "
                            f"do not agree with the bond valence analysis."
                        )
            else:
                text_des.append(
                    "Oxidation states from BVA analyzer cannot be determined. "
                    "Thus BVA charge comparison is not conducted."
                )

        elif key == "dos_comparisons":
            comp_types = []
            tani_index = []
            for orb in val:
                if orb.split("_")[-1] in ["s", "p", "d", "f", "summed"]:
                    comp_types.append(orb.split("_")[-1])
                    tani_index.append(str(val[orb]))
            text_des.append(
                "The Tanimoto index from DOS comparisons in the energy range between {}, {} eV "
                "for {} orbitals are: {}.".format(
                    val["e_range"][0],
                    val["e_range"][1],
                    ", ".join(comp_types),
                    ", ".join(tani_index),
                )
            )

    return text_des
