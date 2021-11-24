import argparse
import json

from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.plotter import CohpPlotter

from cohp.analyze import Analysis
from cohp.describe import Description

parser = argparse.ArgumentParser(description='Add some integers.')

parser.add_argument('--description', action="store_true", default=False, help='will perform bonding analysis study')
parser.add_argument('--automaticplot', action="store_true", default=False,
                    help='plots most important interactions automatically')
parser.add_argument('--plotcohp', dest="plotcohp", nargs='+', default=None, type=int, help='plots specific cohps')

parser.add_argument('--integratecohp', action="store_true", help='integrate plots of most important interactions')

parser.add_argument('--POSCAR', dest="poscar", default="POSCAR", type=str, help='path to POSCAR. Default is "POSCAR"')
parser.add_argument('--ylim', dest="ylim", nargs='+', default=None, type=float, help='energy lim for plots')
parser.add_argument('--xlim', dest="xlim", nargs='+', default=None, type=float, help='COHP lim for plots')
parser.add_argument('--charge', default="CHARGE.lobster", type=str,
                    help='path to Charge.lobster. Default is "CHARGE.lobster"')
parser.add_argument('--icohplist', default="ICOHPLIST.lobster", type=str,
                    help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"')
parser.add_argument('--cohpcar', default="COHPCAR.lobster", type=str,
                    help='path to COHPCAR.lobster. Default is "COHPCAR.lobster"')
parser.add_argument('--json', action="store_true",
                    help='will produce a lobsterpy.json with the most important informations')
parser.add_argument('--filename', default="lobsterpy.json", type=str,
                    help='path to ICOHPLIST.lobster. Default is "ICOHPLIST.lobster"')

args = parser.parse_args()

if args.description or args.automaticplot:
    analyse = Analysis(path_to_poscar=args.poscar, path_to_charge=args.charge, path_to_cohpcar=args.cohpcar,
                       path_to_icohplist=args.icohplist)

if args.description:
    describe = Description(analysis_object=analyse)
    describe.write_description()

if args.automaticplot:
    plt = describe.plot_cohps(ylim=args.ylim, xlim=args.xlim, integrated=args.integratecohp)

if args.json:
    analysedict = analyse.condensed_bonding_analysis
    with open("lobsterpy.json", "w") as fd:
        json.dump(analysedict, fd)

if args.plotcohp:

    completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=args.cohpcar, structure_file=args.poscar)

    cp = CohpPlotter()
    # get a nicer plot label

    for label in args.plotcohp:
        cp.add_cohp(label, completecohp.get_cohp_by_label(label=str(label)))
        # check which COHP you are plotting

    x = cp.get_plot(integrated=args.integratecohp)
    x.ylim(args.ylim)
    x.xlim(args.xlim)

    x.show()
