# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""Script to analyze Lobster outputs from the command line."""
from __future__ import annotations

import argparse
import json
import os
from math import log, sqrt
from pathlib import Path

import matplotlib.style
from monty.os.path import zpath
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.io.lobster import Icohplist

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from lobsterpy.featurize.utils import get_file_paths
from lobsterpy.plotting import (
    IcohpDistancePlotter,
    PlainCohpPlotter,
    PlainDosPlotter,
    get_style_list,
)


def main() -> None:
    """Entry point for setup.py installer."""
    args = get_parser().parse_args()
    run(args)


def get_parser() -> argparse.ArgumentParser:
    """Construct argumentparser with subcommands and sections."""
    parser = argparse.ArgumentParser(description="Analyze and plot results from Lobster runs.")

    # Arguments that are needed by different actions, but not always

    incar_file = argparse.ArgumentParser(add_help=False)
    incar_file.add_argument(
        "-fincar",
        "--file-incar",
        default="INCAR",
        dest="incar",
        type=Path,
        help='path to INCAR. Default is "INCAR".',
    )

    charge_file = argparse.ArgumentParser(add_help=False)
    charge_file.add_argument(
        "-fcharge",
        "--file-charge",
        default="CHARGE.lobster",
        dest="charge",
        type=Path,
        help='path to CHARGE.lobster. Default is "CHARGE.lobster"',
    )

    structure_file = argparse.ArgumentParser(add_help=False)
    structure_file.add_argument(
        "-fstruct",
        "--file-structure",
        default="POSCAR",
        dest="structure",
        type=Path,
        help='path to structure file. Default is "POSCAR". '
        'Can also read "POSCAR.lobster" file or any '
        'suitable file format supported by pymatgen "Structure.from_file" method.',
    )

    potcar_file = argparse.ArgumentParser(add_help=False)
    potcar_file.add_argument(
        "-fpotcar",
        "--file-potcar",
        default="POTCAR",
        dest="potcar",
        type=Path,
        help='path to POTCAR. Default is "POTCAR".',
    )

    icohplist_file = argparse.ArgumentParser(add_help=False)
    icohplist_file.add_argument(
        "-ficohp",
        "--file-icohplist",
        default="ICOHPLIST.lobster",
        dest="icohplist",
        type=Path,
        help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"',
    )

    coxxcar_file = argparse.ArgumentParser(add_help=False)
    coxxcar_file.add_argument(
        "-fcohp",
        "--file-cohpcar",
        default="COHPCAR.lobster",
        dest="cohpcar",
        type=Path,
        help='path to COHPCAR.lobster. Default is "COHPCAR.lobster". This argument '
        "can also read COBICARs or COOPCARs. One needs to use appropriate --cobis or "
        "--coops options along with this argument when plotting",
    )

    doscar_file = argparse.ArgumentParser(add_help=False)
    doscar_file.add_argument(
        "-fdos",
        "--file-doscar",
        default="DOSCAR.lobster",
        dest="doscar",
        type=Path,
        help='path to DOSCAR.lobster. Default is "DOSCAR.lobster".',
    )

    user_basis_arg = argparse.ArgumentParser(add_help=False)
    user_basis_arg.add_argument(
        "--userbasis",
        "--user-basis",
        default=None,
        type=_element_basis,
        nargs="+",
        help="Use the specific basis provided by the user to "
        "generate the inputs (e.g.,  --userbasis Cr.3d.3p.4s N.2s.2p). "
        "Default is None.",
    )

    # groups of arguments for specific actions

    calc_quality_description_file_parent = argparse.ArgumentParser(add_help=False)

    # Input args for specifically needed for calculation quality description
    calc_quality_description_file_group = calc_quality_description_file_parent.add_argument_group("")
    calc_quality_description_file_group.add_argument(
        "-fvasprun",
        "--file-vasprun",
        default="vasprun.xml",
        dest="vasprun",
        type=Path,
        help='path to vasprun.xml. Default is "vasprun.xml".',
    )
    calc_quality_description_file_group.add_argument(
        "-fbandoverlaps",
        "--file-bandoverlaps",
        default="bandOverlaps.lobster",
        dest="bandoverlaps",
        type=Path,
        help='path to bandOverlaps.lobster. Default is "bandOverlaps.lobster".',
    )
    calc_quality_description_file_group.add_argument(
        "-flobsterin",
        "--file-lobsterin",
        default="lobsterin",
        dest="lobsterin",
        type=Path,
        help='path to lobsterin. Default is "lobsterin".',
    )
    calc_quality_description_file_group.add_argument(
        "-flobsterout",
        "--file-lobsterout",
        default="lobsterout",
        dest="lobsterout",
        type=Path,
        help='path to lobsterout. Default is "lobsterout".',
    )
    calc_quality_description_file_group.add_argument(
        "-potsymbols",
        "--potcar-symbols",
        dest="potcarsymbols",
        type=_potcar_symbols,
        # nargs="+",
        help="List of potcar symbols",
    )

    # group of arguments related to writing LOBSTER calcs inputs
    output_parent = argparse.ArgumentParser(add_help=False)
    output_file_group = output_parent.add_argument_group("Output files")
    output_file_group.add_argument(
        "-fincarout",
        "--file-incar-out",
        default="INCAR.lobsterpy",
        dest="incarout",
        type=Path,
        help='path to INCAR that lobsterpy generates. Default is "INCAR.lobsterpy"',
    )
    output_file_group.add_argument(
        "-flobsterin",
        "--file-lobsterin",
        default="lobsterin.lobsterpy",
        dest="lobsterinout",
        type=Path,
        help='Base for path to lobsterins that lobsterpy generates. Default is "lobsterin.lobsterpy"',
    )
    output_file_group.add_argument(
        "--overwrite",
        "--overwrite-files",
        dest="overwrite",
        default=False,
        action="store_true",
        help="overwrites already created INCARs and lobsterins with the give name.",
    )
    # TODO: Add some output arguments: options to supply your own basis
    # General matplotlib plotting arguments common to all kinds of plots

    plotting_parent = argparse.ArgumentParser(add_help=False)
    plotting_group = plotting_parent.add_argument_group("Plotting")
    plotting_group.add_argument(
        "--fontsize",
        "--font-size",
        type=float,
        dest="fontsize",
        default=None,
        help="Base font size",
    )
    plotting_group.add_argument(
        "-ht",
        "--height",
        type=float,
        dest="height",
        default=None,
        help="Plot height in inches",
    )
    plotting_group.add_argument(
        "--hideplot",
        "--hide-plot",
        dest="hideplot",
        action="store_true",
        help="Hide plot output. Especially relevant when plots are saved.",
    )
    plotting_group.add_argument(
        "-nbs",
        "--no-base-style",
        "--nobasestyle",
        action="store_true",
        dest="no_base_style",
        help=(
            "Disable inbuilt style entirely. This may prevent interference with external "
            "stylesheets when using --style."
        ),
    )
    plotting_group.add_argument(
        "-sty",
        "--style",
        type=str,
        nargs="+",
        dest="style",
        default=None,
        help="Matplotlib style sheet(s) for plot appearance",
    )
    plotting_group.add_argument(
        "--save-plot",
        "--saveplot",
        "-s",
        type=Path,
        metavar="FILENAME",
        default=None,
        dest="save_plot",
        help="Save plot to file",
    )
    plotting_group.add_argument("--title", type=str, default="", help="Plot title")
    plotting_group.add_argument(
        "-wd",
        "--width",
        type=float,
        dest="width",
        default=None,
        help="Plot width in inches",
    )
    plotting_group.add_argument(
        "-x",
        "--xlim",
        dest="xlim",
        nargs=2,
        default=None,
        type=float,
        help="Set x-axis limits for the plots",
    )
    plotting_group.add_argument(
        "-y",
        "--ylim",
        dest="ylim",
        nargs=2,
        default=None,
        type=float,
        help="Set y-axis limits for the plots",
    )
    broadening_group = plotting_group.add_mutually_exclusive_group()
    broadening_group.add_argument_group("Broadening")
    broadening_group.add_argument(
        "--fwhm",
        type=float,
        default=None,
        help="Full-width-half-maximum of Gaussian broadening.",
    )
    broadening_group.add_argument(
        "--sigma",
        type=float,
        default=None,
        help="Standard deviation of Gaussian broadening.",
    )
    # Group pf arguments specific to dos plotter
    dos_plotting_parent = argparse.ArgumentParser(add_help=False)
    dos_plotting_group = dos_plotting_parent.add_argument_group("Plotting")
    dos_plotting_group.add_argument(
        "-addtdos",
        "--addtotal",
        "--add-total",
        action="store_true",
        help="Add total dos to DOS plot.",
    )
    dos_plotting_group.add_argument(
        "-el",
        "--element",
        type=str,
        dest="element",
        nargs="+",
        default=None,
        help="Add spd DOS projections for requested element to the DOS plot",
    )
    dos_plotting_group.add_argument(
        "-orb",
        "--orbital",
        type=str,
        dest="orbital",
        nargs="+",
        default=None,
        help="Orbital name for the site for which DOS are to be added",
    )
    dos_plotting_group.add_argument(
        "--site",
        type=int,
        dest="site",
        nargs="+",
        default=None,
        help="Site index in the crystal structure for which DOS need to be added",
    )
    dos_plotting_group.add_argument(
        "-sspin",
        "--summedspins",
        "--summed-spins",
        dest="summedspins",
        action="store_true",
        help="Plot summed spins DOS",
    )

    group = dos_plotting_parent.add_mutually_exclusive_group()
    group.add_argument(
        "-eldos",
        "--elementdos",
        dest="elementdos",
        action="store_true",
        help="Add DOS projections for each element to the DOS plot",
    )
    group.add_argument(
        "--spddos",
        "--spd-dos",
        dest="spddos",
        action="store_true",
        help="Add spd projected dos to the DOS plot",
    )
    # Argument common to COHPs / COOPs / COBIs or DOS plots
    advanced_plotting_args = argparse.ArgumentParser(add_help=False)
    advanced_plotting_args.add_argument(
        "--invertaxis",
        "--invert-axis",
        dest="invertaxis",
        action="store_true",
        help="Invert plot axis of DOS or COOPs COHPs or COBIS",
    )
    # Argument specific to COHPs / COOPs / COBIs plots
    coxx_plotting_args = argparse.ArgumentParser(add_help=False)
    coxx_plotting_args.add_argument(
        "-integr",
        "--integrated",
        dest="integrated",
        action="store_true",
        help="Show integrated cohp/cobi/coop plots.",
    )
    # Arguments specific to lobsterpy.cohp.analyze.Analysis and
    # lobsterpy.cohp.describe.Description class
    auto_parent = argparse.ArgumentParser(add_help=False)
    auto_group = auto_parent.add_argument_group("Adjustable automatic analysis parameters")
    auto_group.add_argument(
        "-allb",
        "--allbonds",
        "--all-bonds",
        action="store_true",
        default=False,
        help="Consider all bonds during the automatic analysis, not only cation-anion bonds (default) ",
    )
    auto_group.add_argument(
        "--cutofficohp",
        type=float,
        default=0.1,
        help="Consider only bonds that are stronger than cutoff_icoxx * strongest ICOXX "
        " (ICOHP or ICOBI or ICOOP) for automatic analysis.",
    )
    auto_group.add_argument(
        "-fjson",
        "--file-json",
        nargs="?",
        type=Path,
        default=None,
        metavar="FILENAME",
        const=Path("lobsterpy.json"),
        help="Write a JSON file with the most important information",
    )
    auto_group.add_argument(
        "--noisecutoff",
        type=float,
        default=None,
        help="Sets the lower limit of icohps or icoops or icobis considered in automatic analysis",
    )
    auto_group.add_argument(
        "-orbresol",
        "--orbitalresolved",
        "--orbital-resolved",
        action="store_true",
        default=False,
        help="Switch on orbital resolved analysis of (I)COHPs or (I)COBIs or (I)COOPs with all relevant orbitals.",
    )
    auto_group.add_argument(
        "-orbcutoff",
        "--orbitalcutoff",
        "--orbital-cutoff",
        type=float,
        default=0.05,
        help="Consider only orbital interactions that are stronger than orbitalintcutoff * 100 of "
        "relevant bonds (ICOHP or ICOBI or ICOOP). Can only be used in combination "
        "with orbital-resolved automatic analysis (--orbitalresolved).",
    )

    auto_group.add_argument(
        "-sspins",
        "--summedspins",
        "--summed-spins",
        action="store_true",
        default=False,
        help="Sum COHP `Spin.up` and `Spin.down` populations for automatic analysis",
    )

    # Argument that will help to switch automatic analysis
    analysis_switch = argparse.ArgumentParser(add_help=False)
    analysis_group = analysis_switch.add_argument_group(
        "Arguments to switch type of files analyzed during automatic analysis"
        " (Also indicates file type for 'plot/plot-icohp-distance' action in cli)"
    )
    analysis_group.add_argument(
        "--cobis",
        "--cobis",
        action="store_true",
        help="Setting this option starts automatic bonding analysis using COBIs"
        " (Also indicates plotter that input file contains COBI data)",
    )
    analysis_group.add_argument(
        "--coops",
        "--coops",
        action="store_true",
        help="Setting this option starts automatic bonding analysis using COOPs"
        " (Also indicates plotter that input file contains COOP data)",
    )
    # Specific to interactive plotter args
    interactive_plotter_args = argparse.ArgumentParser(add_help=False)
    interactive_plotter_group = interactive_plotter_args.add_argument_group("Options specific to interactive plotter")
    interactive_plotter_group.add_argument(
        "-labresol",
        "--labelresolved",
        "--label-resolved",
        action="store_true",
        help="Create automatic interactive plots with all relevant bond labels. "
        "If not set, plots consist of summed cohps.",
    )

    # Args specific to calc quality description dict and texts
    calc_quality_args = argparse.ArgumentParser(add_help=False)
    calc_quality_args_group = calc_quality_args.add_argument_group(
        "Options to change default output of quality description"
    )
    calc_quality_args_group.add_argument(
        "--bvacomp",
        "--bva-comp",
        action="store_true",
        default=False,
        help="Enable the BVA charge comparison in automatic LOBSTER calc quality analysis ",
    )
    calc_quality_args_group.add_argument(
        "--doscomp",
        "--dos-comp",
        action="store_true",
        default=False,
        help="Enable the DOS comparison in automatic LOBSTER calc quality analysis ",
    )
    calc_quality_args_group.add_argument(
        "--erange",
        dest="erange",
        nargs=2,
        default=[-5, 0],
        type=int,
        help="Energy range for DOS comparisons",
    )
    calc_quality_args_group.add_argument(
        "-fcalcqualjson",
        "--file-calc-quality-json",
        nargs="?",
        type=Path,
        default=None,
        metavar="FILENAME",
        const=Path("calc_quality_json.json"),
        help="Write a JSON file with the LOBSTER calc quality analysis",
    )
    calc_quality_args_group.add_argument(
        "--nbins",
        dest="nbins",
        default=None,
        type=int,
        help="Number of bins for DOS comparisons",
    )
    # Build the actions using arguments defined earlier
    subparsers = parser.add_subparsers(
        dest="action",
        required=True,
        help="Use -h/--help after the chosen subcommand to see further options.",
    )
    subparsers.add_parser(
        "create-inputs",
        aliases=["createinputs"],
        parents=[
            incar_file,
            structure_file,
            potcar_file,
            user_basis_arg,
            output_parent,
        ],
        help="Create inputs for lobster computation. It works only with PBE POTCARs, "
        "as, currently, only the pbeVASPfit2015 basis fitted to PBE POTCARs in LOBSTER includes "
        "additional orbitals relevant for the solid-state materials. Please check out our "
        "publication https://doi.org/10.1002/cplu.202200123 and LOBSTER program manual for more information.",
    )
    subparsers.add_parser(
        "description",
        parents=[
            charge_file,
            coxxcar_file,
            icohplist_file,
            structure_file,
            auto_parent,
            analysis_switch,
        ],
        help=("Deliver a text description from automatic analysis of COHPs or COBIS or COOP results from Lobster run"),
    )
    subparsers.add_parser(
        "description-quality",
        parents=[
            calc_quality_description_file_parent,
            charge_file,
            doscar_file,
            potcar_file,
            structure_file,
            calc_quality_args,
        ],
        help=(
            "Deliver a text description of the LOBSTER calc quality analysis. "
            "Mandatory required files: POSCAR, POTCAR, lobsterout, lobsterin. "
            "Optional files (BVA comparison): CHARGE.lobster, "
            "(DOS comparison): DOSCAR.lobster/ DOSCAR.LSO.lobster, Vasprun.xml."
        ),
    )
    subparsers.add_parser(
        "plot-automatic",
        aliases=[
            "plot-auto",
            "automatic-plot",
            "automaticplot",
            "auto-plot",
            "autoplot",
        ],
        parents=[
            charge_file,
            coxxcar_file,
            icohplist_file,
            structure_file,
            auto_parent,
            plotting_parent,
            coxx_plotting_args,
            analysis_switch,
        ],
        help=(
            "Plot most important COHPs or COBIs or COOPs automatically."
            " This option also includes an automatic description."
        ),
    )
    subparsers.add_parser(
        "plot-automatic-ia",
        aliases=[
            "plot-auto-ia",
            "automatic-plot-ia",
            "automaticplotia",
            "autoplotia",
            "auto-plot-ia",
        ],
        parents=[
            charge_file,
            coxxcar_file,
            icohplist_file,
            structure_file,
            auto_parent,
            plotting_parent,
            coxx_plotting_args,
            interactive_plotter_args,
            analysis_switch,
        ],
        help=("Creates an interactive plot of most important COHPs or COBIs or COOPs automatically."),
    )
    subparsers.add_parser(
        "plot-dos",
        aliases=["plotdos"],
        parents=[
            doscar_file,
            structure_file,
            dos_plotting_parent,
            advanced_plotting_args,
            plotting_parent,
        ],
        help="Plots DOS from lobster computation.",
    )
    subparsers.add_parser(
        "plot-icohp-distance",
        aliases=["ploticohpdistance"],
        parents=[icohplist_file, plotting_parent, analysis_switch],
        help="Plot ICOHPs or ICOOPs or ICOBIs with respect to bond lengths",
    )
    # Mode for normal plotting (without automatic detection of relevant COHPs)
    plot_parser = subparsers.add_parser(
        "plot",
        parents=[
            coxxcar_file,
            structure_file,
            plotting_parent,
            advanced_plotting_args,
            coxx_plotting_args,
            analysis_switch,
        ],
        help="Plot specific COHPs/COBIs/COOPs based on bond numbers.",
    )
    plot_parser.add_argument(
        "bond_numbers",
        nargs="+",
        type=int,
        help="List of bond numbers, determining COHPs/COBIs/COOPs to include in plot.",
    )
    plot_grouping = plot_parser.add_mutually_exclusive_group()
    plot_grouping.add_argument(
        "--summed",
        action="store_true",
        help="Show a summed COHP",
    )
    plot_grouping.add_argument(
        "-orbwise",
        "--orbitalwise",
        dest="orbitalwise",
        nargs="+",
        default=None,
        type=str,
        help=(
            "Plot cohps of specific orbitals. e.g. to plot 2s-2s interaction of "
            'bond with label 1, use "lobsterpy plot 1 --orbitalwise 2s-2s". '
            'To plot all orbitalwise cohps of one bond, you can use "all" instead of "2s-2s". '
            "To plot orbitalwise interactions of more than one bond, use, for example, "
            '"lobsterpy plot 1 1 --orbitalwise "3s-3s" "2px-3s"'
        ),
    )
    return parser


def _element_basis(string: str):
    """
    Parse element and basis from string.

    :param string: string to parse

    Returns:
            element, basis
    """
    cut_list = string.split(".")
    element = cut_list[0]
    basis = " ".join(cut_list[1:])
    return element, basis


def _potcar_symbols(string: str):
    """
    Parse string of potcar symbols and return a list.

    :param string: string of potcar symbols

    Returns:
        list of potcar symbols
    """
    return string.split(" ")


def _user_figsize(width, height, aspect=None):
    """Get figsize options from user input, if any.

    If only width x or height is provided, use a target aspect ratio to derive
    the other one.

    Returns a dict which can be merged into style kwargs
    """
    if width is None and height is None:
        return {}
    if width is not None and height is not None:
        return {"figure.figsize": (width, height)}

    if aspect is None:
        aspect = (sqrt(5) + 1) / 2  # Golden ratio
    if width is None:
        return {"figure.figsize": (height * aspect, height)}
    return {"figure.figsize": (width, width / aspect)}


# TODO: add automatic functionality for COBIs, COOPs
def run(args):
    """
    Run actions based on args.

    :param args: args for cli

    """
    if args.action in [
        "plot-auto",
        "automatic-plot",
        "automaticplot",
        "auto-plot",
        "autoplot",
    ]:
        args.action = "plot-automatic"

    if args.action in [
        "plot-auto-ia",
        "automatic-plot-ia",
        "automaticplotia",
        "autoplotia",
        "auto-plot-ia",
    ]:
        args.action = "plot-automatic-ia"

    if args.action in [
        "description",
        "plot-automatic",
        "plot-automatic-ia",
    ]:
        req_files = get_file_paths(
            path_to_lobster_calc=Path(os.getcwd()), requested_files=["structure", "charge", "icohplist", "cohpcar"]
        )

        if args.coops:
            req_files_coops = get_file_paths(
                path_to_lobster_calc=Path(os.getcwd()), requested_files=["icooplist", "coopcar"]
            )

            req_files["icohplist"] = req_files_coops["icooplist"]
            req_files["cohpcar"] = req_files_coops["coopcar"]

        if args.cobis:
            req_files_cobis = get_file_paths(
                path_to_lobster_calc=Path(os.getcwd()), requested_files=["icobilist", "cobicar"]
            )

            req_files["icohplist"] = req_files_cobis["icobilist"]
            req_files["cohpcar"] = req_files_cobis["cobicar"]

        for arg_name in req_files:
            setattr(args, arg_name, req_files[arg_name])

    if args.action in ["description", "plot-automatic", "plot-automatic-ia"]:
        which_bonds = "all" if args.allbonds else "cation-anion"

        analyse = Analysis(
            path_to_poscar=args.structure,
            path_to_charge=args.charge,
            path_to_cohpcar=args.cohpcar,
            path_to_icohplist=args.icohplist,
            which_bonds=which_bonds,
            are_coops=args.coops,
            are_cobis=args.cobis,
            noise_cutoff=args.noisecutoff,
            cutoff_icohp=args.cutofficohp,
            orbital_cutoff=args.orbitalcutoff,
            orbital_resolved=args.orbitalresolved,
            summed_spins=args.summedspins,
        )

        describe = Description(analysis_object=analyse)
        describe.write_description()

        if args.file_json is not None:
            analysedict = analyse.condensed_bonding_analysis
            with open(args.file_json, "w") as fd:
                json.dump(analysedict, fd)

    if args.action in [
        "plot",
        "plot-automatic",
        "plot-automatic-ia",
        "plot-dos",
        "plotdos",
        "plot-icohp-distance",
        "ploticohpdistance",
    ]:
        style_kwargs = {}
        style_kwargs.update(_user_figsize(args.width, args.height))
        if args.fontsize:
            style_kwargs.update({"font.size": args.fontsize})

        style_list = get_style_list(no_base_style=args.no_base_style, styles=args.style, **style_kwargs)
        matplotlib.style.use(style_list)

        if args.sigma:
            sigma = args.sigma
        elif args.fwhm:
            sigma = args.fwhm / (2 * sqrt(2 * log(2)))
        else:
            sigma = None

    if args.action in ["plot-automatic"]:
        describe.plot_cohps(
            ylim=args.ylim,
            xlim=args.xlim,
            integrated=args.integrated,
            save=args.save_plot is not None,
            filename=args.save_plot,
            title=args.title,
            sigma=sigma,
            hide=args.hideplot,
        )

    if args.action in ["plot-automatic-ia"]:
        describe.plot_interactive_cohps(
            ylim=args.ylim,
            xlim=args.xlim,
            integrated=args.integrated,
            save_as_html=args.save_plot,
            filename=args.save_plot,
            title=args.title,
            sigma=sigma,
            hide=args.hideplot,
            label_resolved=args.labelresolved,
            orbital_resolved=args.orbitalresolved,
        )

    if args.action == "plot":
        if args.cobis:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["cobicar"]).get(
                "cobicar"
            )
            options = {"are_cobis": True, "are_coops": False}
        elif args.coops:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["coopcar"]).get(
                "coopcar"
            )
            options = {"are_cobis": False, "are_coops": True}
        else:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["cohpcar"]).get(
                "cohpcar"
            )
            options = {"are_cobis": False, "are_coops": False}

        struture_filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["structure"]).get(
            "structure"
        )

        completecohp = CompleteCohp.from_file(
            fmt="LOBSTER",
            filename=filename,
            structure_file=struture_filename,
            **options,
        )
        cp = PlainCohpPlotter(**options)

        if not args.summed:
            # TODO: add checks for label in allowed labels -> print all labels
            # TODO: add check if args.oribtalwise is exactly as long as labels
            # TODO: add check if orbital is in args.orbitalwise

            for label in args.bond_numbers:
                if str(label) not in completecohp.bonds:
                    raise IndexError(
                        "The provided bond label " + str(label) + " is not available in ICO**LIST.lobster.\n "
                        "Allowed options are in this list: \n"
                        + str([int(listi) for listi in list(completecohp.bonds.keys())])
                    )

            if not args.orbitalwise:
                for label in args.bond_numbers:
                    cp.add_cohp(label, completecohp.get_cohp_by_label(label=str(label)))
            else:
                if len(args.bond_numbers) != len(args.orbitalwise):
                    raise IndexError(
                        "Please provide as mainy orbitals as bond labels,"
                        " e.g., lobsterpy plot 1 1 --orbitalwise '2s-2s' '2s-2px'"
                    )

                for ilabel, label in enumerate(args.bond_numbers):
                    orbitals = args.orbitalwise[ilabel]

                    availableorbitals = list(completecohp.orb_res_cohp[str(label)].keys())
                    orbitaloptions = [*availableorbitals, "all"]

                    if orbitals not in orbitaloptions:
                        raise IndexError(
                            "Orbital in not available for current bond. \n"
                            "For bond "
                            + str(label)
                            + " only the following orbital options are available: \n"
                            + str(orbitaloptions)
                        )

                    if orbitals != "all":
                        cp.add_cohp(
                            str(label) + ": " + orbitals,
                            completecohp.get_orbital_resolved_cohp(label=str(label), orbitals=orbitals),
                        )
                    else:
                        for orbitals in availableorbitals:
                            cp.add_cohp(
                                str(label) + ": " + orbitals,
                                completecohp.get_orbital_resolved_cohp(label=str(label), orbitals=orbitals),
                            )
        else:
            cp.add_cohp(
                str(args.bond_numbers),
                completecohp.get_summed_cohp_by_label_list(label_list=[str(label) for label in args.bond_numbers]),
            )

        plt = cp.get_plot(integrated=args.integrated, xlim=args.xlim, ylim=args.ylim, sigma=sigma)

        ax = plt.gca()
        ax.set_title(args.title)

        if not args.hideplot and not args.save_plot:
            plt.show()
        elif args.save_plot and not args.hideplot:
            plt.show()
            fig = plt.gcf()
            fig.savefig(args.save_plot)
        if args.save_plot and args.hideplot:
            plt.savefig(args.save_plot)

    if args.action in ["create-inputs", "createinputs"]:
        from pymatgen.core.structure import Structure
        from pymatgen.io.lobster import Lobsterin

        # Check for .gz files exist for default values and update accordingly
        default_files = {
            "structure": "POSCAR",
            "potcar": "POTCAR",
            "incar": "INCAR",
        }

        for arg_name in default_files:
            file_path = getattr(args, arg_name)
            if not file_path.exists():
                gz_file_path = file_path.with_name(zpath(file_path.name))
                if gz_file_path.exists():
                    setattr(args, arg_name, gz_file_path)
                else:
                    raise ValueError(
                        "Files necessary for creating inputs for LOBSTER calcs not found in the current directory."
                    )

        if args.userbasis is None:
            # This will rely on standard basis files as stored in pymatgen

            potcar_names = Lobsterin._get_potcar_symbols(POTCAR_input=args.potcar)

            list_basis_dict = Lobsterin.get_all_possible_basis_functions(
                structure=Structure.from_file(args.structure),
                potcar_symbols=potcar_names,
            )

            for ibasis, basis_dict in enumerate(list_basis_dict):
                lobsterinput = Lobsterin.standard_calculations_from_vasp_files(
                    args.structure,
                    args.incar,
                    None,
                    option="standard",
                    dict_for_basis=basis_dict,
                )

                lobsterin_path = Path(str(args.lobsterinout) + "-" + str(ibasis))
                incar_path = Path(str(args.incarout) + "-" + str(ibasis))

                if (not lobsterin_path.is_file() and not incar_path.is_file()) or (args.overwrite):
                    lobsterinput.write_lobsterin(lobsterin_path)
                    lobsterinput.write_INCAR(
                        incar_input=args.incar,
                        incar_output=incar_path,
                        poscar_input=args.structure,
                        isym=0,
                    )
                else:
                    raise ValueError('please use "--overwrite" if you would like to overwrite existing lobster inputs')
        else:
            # convert list userbasis to dict
            userbasis = {}
            for userbasis_single in args.userbasis:
                userbasis[userbasis_single[0]] = userbasis_single[1]

            lobsterinput = Lobsterin.standard_calculations_from_vasp_files(
                args.structure,
                args.incar,
                None,
                option="standard",
                dict_for_basis=userbasis,
            )

            lobsterin_path = Path(str(args.lobsterinout) + "-" + str(0))
            incar_path = Path(str(args.incarout) + "-" + str(0))

            if (not lobsterin_path.is_file() and not incar_path.is_file()) or (args.overwrite):
                lobsterinput.write_lobsterin(lobsterin_path)
                lobsterinput.write_INCAR(
                    incar_input=args.incar,
                    incar_output=incar_path,
                    poscar_input=args.structure,
                    isym=0,
                )
            else:
                raise ValueError('please use "--overwrite" if you would like to overwrite existing lobster inputs')

    if args.action in ["description-quality"]:
        # Check for .gz files exist for default values and update accordingly
        req_files = get_file_paths(
            path_to_lobster_calc=Path(os.getcwd()), requested_files=["structure", "lobsterin", "lobsterout"]
        )
        for arg_name in req_files:
            setattr(args, arg_name, req_files[arg_name])

        optional_files = {
            "bandoverlaps": "bandOverlaps.lobster",
            "potcar": "POTCAR",
        }

        for arg_name in optional_files:
            file_path = getattr(args, arg_name)
            if not file_path.exists():
                gz_file_path = file_path.with_name(zpath(file_path.name))
                if gz_file_path.exists():
                    setattr(args, arg_name, gz_file_path)

        bva_comp = args.bvacomp

        if bva_comp:
            bva_files = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["charge"])
            for arg_name in bva_files:
                setattr(args, arg_name, bva_files[arg_name])

        dos_comparison = args.doscomp

        if dos_comparison:
            dos_files = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["doscar", "vasprun"])
            for arg_name in dos_files:
                if arg_name == "doscar":
                    file_path = getattr(args, arg_name)
                    if not file_path.exists():
                        gz_file_path = file_path.with_name(zpath(file_path.name))
                        if gz_file_path.exists():
                            setattr(args, arg_name, gz_file_path)
                        else:
                            raise ValueError(f"{file_path} or {gz_file_path} not found in current directory")
                else:
                    setattr(args, arg_name, dos_files[arg_name])

        potcar_file_path = args.potcar

        quality_dict = Analysis.get_lobster_calc_quality_summary(
            path_to_poscar=args.structure,
            path_to_charge=args.charge,
            path_to_lobsterout=args.lobsterout,
            path_to_lobsterin=args.lobsterin,
            path_to_potcar=None if not potcar_file_path.exists() else potcar_file_path,
            potcar_symbols=args.potcarsymbols,
            path_to_bandoverlaps=args.bandoverlaps,
            dos_comparison=dos_comparison,
            bva_comp=bva_comp,
            path_to_doscar=args.doscar,
            e_range=args.erange,
            n_bins=args.nbins,
            path_to_vasprun=args.vasprun,
        )

        quality_text = Description.get_calc_quality_description(quality_dict)
        Description.write_calc_quality_description(quality_text)

        if args.file_calc_quality_json is not None:
            with open(args.file_calc_quality_json, "w") as fd:
                json.dump(quality_dict, fd)

    if args.action in ["plot-dos", "plotdos"]:
        req_files = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["structure", "doscar"])

        for arg_name in req_files:
            if arg_name == "doscar":
                file_path = getattr(args, arg_name)
                if not file_path.exists():
                    gz_file_path = file_path.with_name(zpath(file_path.name))
                    if gz_file_path.exists():
                        setattr(args, arg_name, gz_file_path)
                    else:
                        raise ValueError(f"{file_path} or {gz_file_path} not found in current directory")
            else:
                setattr(args, arg_name, req_files[arg_name])

        from pymatgen.io.lobster import Doscar

        lobs_dos = Doscar(doscar=args.doscar, structure_file=args.structure).completedos

        dos_plotter = PlainDosPlotter(summed=args.summedspins, sigma=args.sigma)
        if args.addtotal:
            dos_plotter.add_dos(dos=lobs_dos, label="Total DOS")
        if args.spddos:
            dos_plotter.add_dos_dict(dos_dict=lobs_dos.get_spd_dos())

        if args.elementdos:
            dos_plotter.add_dos_dict(dos_dict=lobs_dos.get_element_dos())

        if args.element:
            for element in args.element:
                element_spddos = lobs_dos.get_element_spd_dos(el=element)
                for orbital, dos in element_spddos.items():
                    label = f"{element}: {orbital.name}"
                    dos_plotter.add_dos_dict(dos_dict={label: dos})

        if args.site is not None and args.orbital:
            if len(args.site) > len(args.orbital):
                for site in args.site:
                    for orbital in args.orbital:
                        dos_plotter.add_site_orbital_dos(site_index=site, orbital=orbital, dos=lobs_dos)
            elif len(args.orbital) > len(args.site):
                for orbital in args.orbital:
                    for site in args.site:
                        dos_plotter.add_site_orbital_dos(site_index=site, orbital=orbital, dos=lobs_dos)
            else:
                for site, orbital in zip(args.site, args.orbital):
                    dos_plotter.add_site_orbital_dos(site_index=site, orbital=orbital, dos=lobs_dos)
        elif (
            args.site is None
            and not args.orbital
            and not args.element
            and not args.spddos
            and not args.elementdos
            and not args.addtotal
        ):
            dos_plotter.add_dos(dos=lobs_dos, label="Total DOS")
            dos_plotter.add_dos_dict(dos_dict=lobs_dos.get_element_dos())

        elif (args.site is None or not args.orbital) and (not args.element and not args.spddos and not args.elementdos):
            raise ValueError("Please set both args i.e site and orbital to generate the plot")

        plt = dos_plotter.get_plot(
            xlim=args.xlim,
            ylim=args.ylim,
            beta_dashed=True,
            invert_axes=args.invertaxis,
        )

        ax = plt.gca()
        ax.set_title(args.title)

        if not args.hideplot and not args.save_plot:
            plt.show()
        elif args.save_plot and not args.hideplot:
            plt.show()
            fig = plt.gcf()
            fig.savefig(args.save_plot)
        if args.save_plot and args.hideplot:
            plt.savefig(args.save_plot)

    if args.action in ["plot-icohp-distance", "ploticohpdistance"]:
        if args.cobis:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["icobilist"]).get(
                "icobilist"
            )
            args.icohplist = filename
            # filename = args.icohplist.parent / "ICOBILIST.lobster"
            # if not filename.exists():
            #     filename = filename.with_name(zpath(filename.name))
            options = {"are_cobis": True, "are_coops": False}
        elif args.coops:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["icooplist"]).get(
                "icooplist"
            )
            args.icohplist = filename
            options = {"are_cobis": False, "are_coops": True}
        else:
            filename = get_file_paths(path_to_lobster_calc=Path(os.getcwd()), requested_files=["icohplist"]).get(
                "icohplist"
            )
            args.icohplist = filename
            options = {"are_cobis": False, "are_coops": False}

        icohpcollection = Icohplist(filename=args.icohplist, **options).icohpcollection
        icohp_plotter = IcohpDistancePlotter(**options)

        icohp_plotter.add_icohps(icohpcollection=icohpcollection, label="")

        plt = icohp_plotter.get_plot(xlim=args.xlim, ylim=args.ylim)

        ax = plt.gca()
        ax.set_title(args.title)

        if not args.hideplot and not args.save_plot:
            plt.show()
        elif args.save_plot and not args.hideplot:
            plt.show()
            fig = plt.gcf()
            fig.savefig(args.save_plot)
        if args.save_plot and args.hideplot:
            plt.savefig(args.save_plot)


if __name__ == "__main__":
    main()
