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


parser = argparse.ArgumentParser(
    description="Analyze and plot results from Lobster runs."
)

# Options for Automatic Analysis
parser.add_argument(
    "--description",
    action="store_true",
    default=False,
    help="This will deliver a text description of the COHP results from Lobster and VASP. Implementation of COBIs and COOPs will follow.",
)
parser.add_argument(
    "--automaticplot",
    action="store_true",
    default=False,
    help="This will plot most important COHPs automatically. Implementation of COBIs and COOPs will follow.",
)
parser.add_argument(
    "--json",
    action="store_true",
    help="This will produce a lobsterpy.json with the most important information",
)
parser.add_argument(
    "--filenamejson",
    default="lobsterpy.json",
    type=Path,
    help="Path to json file storing the most important bonding information from the automatic analysis. Default is lobsterpy.json",
)
parser.add_argument(
    "--allbonds",
    action="store_true",
    default=False,
    help="This option will force the automatc analysis to consider all bonds, not only cation-anion bonds (default) ",
)

# options for normal plotting (without automatic detection of relevant COHPs)
parser.add_argument(
    "--plot",
    dest="plot",
    nargs="+",
    default=None,
    type=int,
    help='This plots specific cohps, cobis, coops based on bond numbers, list them after --plot (e.g., "--plot 1"). Default is a COHP plot. You cannot use --plot at the same time as --automaticplot or --description.',
)
parser.add_argument(
    "--cobis",
    "--cobi",
    action="store_true",
    help="if --plot is used as well, it will plot cobis",
)
parser.add_argument(
    "--coops",
    "--coop",
    action="store_true",
    help="if --plot is used as well, it will plot coops",
)
parser.add_argument(
    "--summed",
    action="store_true",
    help='if --plot is used as well, then a summed COHP is shown. Usage: "--plot 1 2 --summed. Cannot be used together with --orbitalwise',
)
parser.add_argument(
    "--orbitalwise",
    dest="orbitalwise",
    nargs="+",
    default=None,
    type=str,
    help='plots cohps of specific orbitals. To plot 2s-2s interaction of bond with label 1, you have to type "lobterpy --plot 1 --orbitalwise 2s-2s". To plot all orbitalwise cohps of one bond, you can use "all" instead of "2s-2s". It cannot be used together with summed at the moment.',
)

# Options for plots
parser.add_argument(
    "--ylim",
    dest="ylim",
    nargs="+",
    default=None,
    type=float,
    help="Energy range for plots",
)
parser.add_argument(
    "--xlim",
    dest="xlim",
    nargs="+",
    default=None,
    type=float,
    help="COHP/COBI/COOP range for plots",
)

# Options for all analysis that can be done with lobsterpy
parser.add_argument(
    "--integrated", action="store_true", help="Show integrated cohp/cobi/coop plots."
)
parser.add_argument(
    "--POSCAR",
    "--poscar",
    dest="poscar",
    default="POSCAR",
    type=Path,
    help='path to POSCAR. Default is "POSCAR"',
)
parser.add_argument(
    "--charge",
    default="CHARGE.lobster",
    type=Path,
    help='path to Charge.lobster. Default is "CHARGE.lobster"',
)
parser.add_argument(
    "--icohplist",
    default="ICOHPLIST.lobster",
    type=Path,
    help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"',
)
parser.add_argument(
    "--cohpcar",
    default="COHPCAR.lobster",
    type=Path,
    help='path to COHPCAR.lobster. Default is "COHPCAR.lobster". This argument will also be read when COBICARs or COOPCARs are plotted.',
)
parser.add_argument(
    "--style",
    type=str,
    nargs="+",
    default=None,
    help="Matplotlib style sheet(s) for plot appearance",
)
parser.add_argument(
    "--no-base-style",
    action="store_true",
    dest="no_base_style",
    help=(
        "Disable inbuilt style entirely. This may prevent interference with external "
        "stylesheets when using --style."
    ),
)
parser.add_argument("--title", type=str, default="", help="Plot title")
parser.add_argument(
    "--save-plot",
    "-s",
    type=str,
    default=None,
    dest="save_plot",
    help="Save plot to file",
)
parser.add_argument("--width", type=float, default=None, help="Plot width in inches")
parser.add_argument("--height", type=float, default=None, help="Plot height in inches")
parser.add_argument("--fontsize", type=float, default=None, help="Base font size")

args = parser.parse_args()


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
def main():
    """
    main of cli
    """
    if args.description or args.automaticplot:
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

    if args.description or args.automaticplot:
        describe = Description(analysis_object=analyse)
        describe.write_description()

    if args.plot or args.automaticplot:
        style_kwargs = {}
        style_kwargs.update(_user_figsize(args.width, args.height))
        if args.fontsize:
            style_kwargs.update({"font.size": args.fontsize})

        style_list = get_style_list(
            no_base_style=args.no_base_style, styles=args.style, **style_kwargs
        )
        matplotlib.style.use(style_list)

    if args.automaticplot:
        describe.plot_cohps(ylim=args.ylim, xlim=args.xlim, integrated=args.integrated)

    if args.json:
        analysedict = analyse.condensed_bonding_analysis
        with open(args.filenamejson, "w") as fd:
            json.dump(analysedict, fd)

    if args.plot and not (args.description or args.automaticplot):
        if args.cobis and args.coops:
            raise ValueError("COBI and COBI cannot be chosen at the same time.")
        if (not args.cobis) and (not args.coops):
            completecohp = CompleteCohp.from_file(
                fmt="LOBSTER", filename=args.cohpcar, structure_file=args.poscar
            )
        else:
            if args.cohpcar.name == "COHPCAR.lobster":
                if args.cobis:
                    filename = args.cohpcar.parent / "COBICAR.lobster"
                elif args.coops:
                    filename = args.cohpcar.parent / "COOPCAR.lobster"
            if args.cobis:
                completecohp = CompleteCohp.from_file(
                    fmt="LOBSTER",
                    filename=filename,
                    structure_file=args.poscar,
                    are_cobis=True,
                )
            elif args.coops:
                completecohp = CompleteCohp.from_file(
                    fmt="LOBSTER",
                    filename=filename,
                    structure_file=args.poscar,
                    are_coops=True,
                )
        if (not args.cobis) and (not args.coops):
            cp = PlainCohpPlotter()
        else:
            if args.cobis:
                cp = PlainCohpPlotter(are_cobis=True)
            elif args.coops:
                cp = PlainCohpPlotter(are_coops=True)

        if not args.summed:
            if not args.orbitalwise:
                for label in args.plot:
                    cp.add_cohp(label, completecohp.get_cohp_by_label(label=str(label)))
            else:
                for ilabel, label in enumerate(args.plot):

                    orbitals = args.orbitalwise[ilabel]
                    if orbitals != "all":
                        cp.add_cohp(
                            str(label) + ": " + orbitals,
                            completecohp.get_orbital_resolved_cohp(
                                label=str(label), orbitals=orbitals
                            ),
                        )
                    else:
                        for orbitals in completecohp.orb_res_cohp[str(label)].keys():
                            cp.add_cohp(
                                str(label) + ": " + orbitals,
                                completecohp.get_orbital_resolved_cohp(
                                    label=str(label), orbitals=orbitals
                                ),
                            )
        else:
            cp.add_cohp(
                str(args.plot),
                completecohp.get_summed_cohp_by_label_list(
                    label_list=[str(label) for label in args.plot]
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
