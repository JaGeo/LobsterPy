# Copyright (c) lobsterpy development team
# Distributed under the terms of a BSD 3-Clause "New" or "Revised" License

"""
Script to analyze Lobster outputs from the command line
"""

import argparse
import json
from math import sqrt
from pathlib import Path

import matplotlib.style

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from lobsterpy.plotting import get_style_list, PlainCohpPlotter
from pymatgen.electronic_structure.cohp import CompleteCohp


def main() -> None:
    """Entry point for setup.py installer"""
    args = get_parser().parse_args()
    run(args)


def get_parser() -> argparse.ArgumentParser:
    """Construct argumentparser with subcommands and sections"""
    parser = argparse.ArgumentParser(
        description="Analyze and plot results from Lobster runs."
    )

    # Arguments common to all actions
    base_parent = argparse.ArgumentParser(add_help=False)
    base_parent.add_argument(
        "--integrated",
        action="store_true",
        help="Show integrated cohp/cobi/coop plots.",
    )
    input_file_group = base_parent.add_argument_group("Input files")
    input_file_group.add_argument(
        "--poscar",
        "--POSCAR",
        dest="poscar",
        default="POSCAR",
        type=Path,
        help='path to POSCAR. Default is "POSCAR"',
    )
    input_file_group.add_argument(
        "--charge",
        default="CHARGE.lobster",
        type=Path,
        help='path to Charge.lobster. Default is "CHARGE.lobster"',
    )
    input_file_group.add_argument(
        "--icohplist",
        default="ICOHPLIST.lobster",
        type=Path,
        help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"',
    )
    input_file_group.add_argument(
        "--cohpcar",
        default="COHPCAR.lobster",
        type=Path,
        help=(
            'path to COHPCAR.lobster. Default is "COHPCAR.lobster". This argument '
            "will also be read when COBICARs or COOPCARs are plotted."
        ),
    )

    plotting_parent = argparse.ArgumentParser(add_help=False)
    plotting_group = plotting_parent.add_argument_group("Plotting")
    plotting_group.add_argument(
        "--ylim",
        dest="ylim",
        nargs=2,
        default=None,
        type=float,
        help="Energy range for plots",
    )
    plotting_group.add_argument(
        "--xlim",
        dest="xlim",
        nargs=2,
        default=None,
        type=float,
        help="COHP/COBI/COOP range for plots",
    )

    plotting_group.add_argument(
        "--style",
        type=str,
        nargs="+",
        default=None,
        help="Matplotlib style sheet(s) for plot appearance",
    )
    plotting_group.add_argument(
        "--no-base-style",
        "--nobasestyle",
        action="store_true",
        dest="no_base_style",
        help=(
            "Disable inbuilt style entirely. This may prevent interference with external "
            "stylesheets when using --style."
        ),
    )
    plotting_group.add_argument("--title", type=str, default="", help="Plot title")
    plotting_group.add_argument(
        "--save-plot",
        "--saveplot",
        "-s",
        type=str,
        metavar="FILENAME",
        default=None,
        dest="save_plot",
        help="Save plot to file",
    )
    plotting_group.add_argument(
        "--width", type=float, default=None, help="Plot width in inches"
    )
    plotting_group.add_argument(
        "--height", type=float, default=None, help="Plot height in inches"
    )
    plotting_group.add_argument(
        "--fontsize", "--font-size", type=float, default=None, help="Base font size"
    )

    auto_parent = argparse.ArgumentParser(add_help=False)
    auto_group = auto_parent.add_argument_group("Automatic analysis")
    auto_group.add_argument(
        "--json",
        nargs="?",
        type=Path,
        default=None,
        metavar="FILENAME",
        const=Path("lobsterpy.json"),
        help="Write a JSON file with the most important information",
    )
    auto_group.add_argument(
        "--allbonds",
        "--all-bonds",
        action="store_true",
        default=False,
        help="This option will force the automatc analysis to consider all bonds, not only cation-anion bonds (default) ",
    )

    subparsers = parser.add_subparsers(
        dest="action",
        required=True,
        help="Use -h/--help after the chosen subcommand to see further options.",
    )
    subparsers.add_parser(
        "description",
        parents=[base_parent, auto_parent],
        help=(
            "Deliver a text description of the COHP results from Lobster "
            "and VASP. Implementation of COBIs and COOPs will follow."
        ),
    )

    subparsers.add_parser(
        "automatic-plot",
        aliases=["automaticplot"],
        parents=[base_parent, auto_parent, plotting_parent],
        help=(
            "Plot most important COHPs automatically. Implementation "
            "of COBIs and COOPs will follow. This option also includes an automatic description."
        ),
    )

    # Mode for normal plotting (without automatic detection of relevant COHPs)
    plot_parser = subparsers.add_parser(
        "plot",
        parents=[base_parent, plotting_parent],
        help="Plot specific COHPs/COBIs/COOPs based on bond numbers.",
    )

    plot_parser.add_argument(
        "bond_numbers",
        nargs="+",
        type=int,
        help="List of bond numbers, determining COHPs/COBIs/COOPs to include in plot.",
    )
    plot_coops_cobis = plot_parser.add_mutually_exclusive_group()
    plot_coops_cobis.add_argument(
        "--cobis",
        "--cobi",
        action="store_true",
        help="Plot COBIs",
    )
    plot_coops_cobis.add_argument(
        "--coops",
        "--coop",
        action="store_true",
        help="Plot COOPs",
    )
    plot_grouping = plot_parser.add_mutually_exclusive_group()
    plot_grouping.add_argument(
        "--summed",
        action="store_true",
        help="Show a summed COHP",
    )
    plot_grouping.add_argument(
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


def _user_figsize(width, height, aspect=None):
    """Get figsize options from user input, if any

    If only width xor height is provided, use a target aspect ratio to derive
    the other one.

    Returns a dict which can be merged into style kwargs
    """

    if width is None and height is None:
        return {}
    elif width is not None and height is not None:
        return {"figure.figsize": (width, height)}
    else:
        if aspect is None:
            aspect = (sqrt(5) + 1) / 2  # Golden ratio
        if width is None:
            return {"figure.figsize": (height * aspect, height)}
        else:
            return {"figure.figsize": (width, width / aspect)}


# TODO: add automatic functionality for COBIs, COOPs
def run(args):
    """

    Args:
        args: args for cli

    Returns:

    """
    if args.action == "automaticplot":
        args.action = "automatic-plot"

    if args.action in ["description", "automatic-plot"]:
        if args.allbonds:
            whichbonds = "all"
        else:
            whichbonds = "cation-anion"
        analyse = Analysis(
            path_to_poscar=args.poscar,
            path_to_charge=args.charge,
            path_to_cohpcar=args.cohpcar,
            path_to_icohplist=args.icohplist,
            whichbonds=whichbonds,
        )

        describe = Description(analysis_object=analyse)
        describe.write_description()

        if args.json is not None:
            analysedict = analyse.condensed_bonding_analysis
            with open(args.json, "w") as fd:
                json.dump(analysedict, fd)

    if args.action in ["plot", "automatic-plot"]:
        style_kwargs = {}
        style_kwargs.update(_user_figsize(args.width, args.height))
        if args.fontsize:
            style_kwargs.update({"font.size": args.fontsize})

        style_list = get_style_list(
            no_base_style=args.no_base_style, styles=args.style, **style_kwargs
        )
        matplotlib.style.use(style_list)

    if args.action in ["automatic-plot"]:
        describe.plot_cohps(
            ylim=args.ylim,
            xlim=args.xlim,
            integrated=args.integrated,
            save=args.save_plot,
            title=args.title,
        )

    if args.action == "plot":
        if args.cobis:
            filename = args.cohpcar.parent / "COBICAR.lobster"
            options = dict(are_cobis=True, are_coops=False)
        elif args.coops:
            filename = args.cohpcar.parent / "COOPCAR.lobster"
            options = dict(are_cobis=False, are_coops=True)
        else:
            filename = args.cohpcar
            options = dict(are_cobis=False, are_coops=False)

        completecohp = CompleteCohp.from_file(
            fmt="LOBSTER", filename=filename, structure_file=args.poscar, **options
        )
        cp = PlainCohpPlotter(**options)

        if not args.summed:
            # TODO: add checks for label in allewod labels -> print all labels
            # TODO: add check if args.oribtalwise is exactly as long as labels
            # TODO: add check if orbital is in args.orbitalwise

            for label in args.bond_numbers:
                if not str(label) in completecohp.bonds.keys():
                    raise IndexError(
                        "The provided bond label "
                        + str(label)
                        + " is not available in ICO**LIST.lobster.\n "
                        "Allowed options are in this list: \n"
                        + str([int(listi) for listi in list(completecohp.bonds.keys())])
                    )

            if not args.orbitalwise:
                for label in args.bond_numbers:
                    cp.add_cohp(label, completecohp.get_cohp_by_label(label=str(label)))
            else:
                if len(args.bond_numbers) != len(args.orbitalwise):
                    raise IndexError(
                        "Please provide as mainy orbitals as bond labels, e.g., lobsterpy plot 1 1 --orbitalwise '2s-2s' '2s-2px'"
                    )

                for ilabel, label in enumerate(args.bond_numbers):
                    orbitals = args.orbitalwise[ilabel]

                    availableorbitals = list(
                        completecohp.orb_res_cohp[str(label)].keys()
                    )
                    orbitaloptions = availableorbitals + ["all"]

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
                            completecohp.get_orbital_resolved_cohp(
                                label=str(label), orbitals=orbitals
                            ),
                        )
                    else:
                        for orbitals in availableorbitals:
                            cp.add_cohp(
                                str(label) + ": " + orbitals,
                                completecohp.get_orbital_resolved_cohp(
                                    label=str(label), orbitals=orbitals
                                ),
                            )
        else:
            cp.add_cohp(
                str(args.bond_numbers),
                completecohp.get_summed_cohp_by_label_list(
                    label_list=[str(label) for label in args.bond_numbers]
                ),
            )

        plt = cp.get_plot(integrated=args.integrated, xlim=args.xlim, ylim=args.ylim)

        ax = plt.gca()
        ax.set_title(args.title)

        if args.save_plot is None:
            plt.show()
        else:
            fig = plt.gcf()
            fig.savefig(args.save_plot)


if __name__ == "__main__":
    main()
