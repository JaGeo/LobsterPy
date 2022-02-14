import argparse
import json
from pathlib import Path

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.plotter import CohpPlotter

parser = argparse.ArgumentParser(description='Analyze and plot results from Lobster runs.')
#Options for Automatic Analysis
parser.add_argument('--description', action="store_true", default=False,
                    help='This will deliver a text description of the COHP results from Lobster and VASP. Implementation of COBIs and COOPs will follow.')
parser.add_argument('--automaticplot', action="store_true", default=False,
                    help='Plots most important COHPs automatically. Implementation of COBIs and COOPs will follow.')
parser.add_argument('--json', action="store_true",
                    help='will produce a lobsterpy.json with the most important information')
parser.add_argument('--filenamejson', default="lobsterpy.json", type=Path,
                    help='path to json file storing the most important bonding information. Default is lobsterpy.json')
parser.add_argument('--allbonds', action="store_true", default=False,
                    help='This option will force the automatc analysis to consider all bonds, not only cation-anion bonds (default) ')

#options for normal plotting (without automatic detection of relevant COHPs)
parser.add_argument('--plot', dest="plot", nargs='+', default=None, type=int,
                    help='Plots specific cohps, cobis, coops based on bond numbers, list them after --plot (e.g., "--plot 1"). Default is a COHP plot. You cannot use --plot at the same time as --automaticplot or --description.')
parser.add_argument('--cobis', '--cobi', action="store_true", help='if --plot is used as well, it will plot cobis')
parser.add_argument('--coops','--coop', action="store_true", help='if --plot is used as well, it will plot coops')
parser.add_argument('--summed', action="store_true",
                    help='if --plot is used as well, then a summed COHP is shown. Usage: "--plot 1 2 --summed. Cannot be used together with --orbitalwise')
parser.add_argument('--orbitalwise', dest="orbitalwise", nargs='+', default=None, type=str,
                    help='plots cohps of specific orbitals. To plot 2s-2s interaction of bond with label 1, you have to type "lobterpy --plot 1 --orbitalwise 2s-2s". To plot all orbitalwise cohps of one bond, you can use "all" instead of "2s-2s". It cannot be used together with summed at the moment.')

#Options for plots
parser.add_argument('--ylim', dest="ylim", nargs='+', default=None, type=float, help='energy limit for plots')
parser.add_argument('--xlim', dest="xlim", nargs='+', default=None, type=float, help='COHP limit for plots')

#Options for all analysis that can be done with lobsterpy
parser.add_argument('--integrated', action="store_true", help='integrate cohp/cobi/coop plots.')
parser.add_argument('--POSCAR', '--poscar', dest="poscar", default="POSCAR", type=Path,
                    help='path to POSCAR. Default is "POSCAR"')
parser.add_argument('--charge', default="CHARGE.lobster", type=Path,
                    help='path to Charge.lobster. Default is "CHARGE.lobster"')
parser.add_argument('--icohplist', default="ICOHPLIST.lobster", type=Path,
                    help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"')
parser.add_argument('--cohpcar', default="COHPCAR.lobster", type=Path,
                    help='path to COHPCAR.lobster. Default is "COHPCAR.lobster". This argument will also be read when COBICARs or COOPCARs are plotted.')


args = parser.parse_args()


# TODO: add automatic functionality for COBIs, COOPs
def main():

    if args.description or args.automaticplot:
        if args.allbonds:
            whichbonds = "all"
        else:
            whichbonds = "cation-anion"
        analyse = Analysis(path_to_poscar=args.poscar, path_to_charge=args.charge, path_to_cohpcar=args.cohpcar,
                           path_to_icohplist=args.icohplist, whichbonds=whichbonds)

    if args.description or args.automaticplot:
        describe = Description(analysis_object=analyse)
        describe.write_description()

    if args.automaticplot:
        plt = describe.plot_cohps(ylim=args.ylim, xlim=args.xlim, integrated=args.integrated)

    if args.json:
        analysedict = analyse.condensed_bonding_analysis
        with open(args.filenamejson, "w") as fd:
            json.dump(analysedict, fd)

    if args.plot and not (args.description or args.automaticplot):
        if args.cobis and args.coops:
            raise ValueError("COBI and COBI cannot be chosen at the same time.")
        if (not args.cobis) and (not args.coops):
            completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=args.cohpcar, structure_file=args.poscar)
        else:
            if args.cohpcar.name == "COHPCAR.lobster":
                if args.cobis:
                    filename = args.cohpcar.parent / "COBICAR.lobster"
                elif args.coops:
                    filename = args.cohpcar.parent / "COOPCAR.lobster"
            if args.cobis:
                completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=filename, structure_file=args.poscar,
                                                      are_cobis=True)
            elif args.coops:
                completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=filename, structure_file=args.poscar,
                                                      are_coops=True)
        if (not args.cobis) and (not args.coops):
            cp = CohpPlotter()
        else:
            if args.cobis:
                cp = CohpPlotter(are_cobis=True)
            elif args.coops:
                cp = CohpPlotter(are_coops=True)

        if not args.summed:
            if not args.orbitalwise:
                for label in args.plot:
                    cp.add_cohp(label, completecohp.get_cohp_by_label(label=str(label)))
            else:
                for ilabel, label in enumerate(args.plot):

                    orbitals = args.orbitalwise[ilabel]
                    if orbitals != "all":
                        cp.add_cohp(str(label) + ': ' + orbitals,
                                    completecohp.get_orbital_resolved_cohp(label=str(label), orbitals=orbitals))
                    else:
                        for orbitals in completecohp.orb_res_cohp[str(label)].keys():
                            cp.add_cohp(str(label) + ': ' + orbitals,
                                        completecohp.get_orbital_resolved_cohp(label=str(label), orbitals=orbitals))
        else:
            cp.add_cohp(str(args.plot),
                        completecohp.get_summed_cohp_by_label_list(label_list=[str(label) for label in args.plot]))

        x = cp.get_plot(integrated=args.integrated)
        x.ylim(args.ylim)
        x.xlim(args.xlim)

        x.show()


if __name__ == "__main__":
    main()
