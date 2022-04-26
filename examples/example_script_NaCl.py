import os

from lobsterpy.cohp.analyze import Analysis
from lobsterpy.cohp.describe import Description

directory = "NaCl"

# Setup analysis dict
analyse = Analysis(
    path_to_poscar=os.path.join(directory, "POSCAR"),
    path_to_icohplist=os.path.join(directory, "ICOHPLIST.lobster"),
    path_to_cohpcar=os.path.join(directory, "COHPCAR.lobster"),
    path_to_charge=os.path.join(directory, "CHARGE.lobster"),
    whichbonds="cation-anion",
)

# Setup Desciption dict
describe = Description(analysis_object=analyse)
describe.write_description()

# Automatic plots
describe.plot_cohps(ylim=[-10, 2], xlim=[-4, 4])

# different dicts that summarize the results
print(analyse.condensed_bonding_analysis)
print(analyse.final_dict_bonds)
print(analyse.final_dict_ions)
