import os
import warnings

from lobsterpy.coxx.analyze import Analysis
from lobsterpy.coxx.describe import Description

warnings.simplefilter("once")

directory = "NaCl"

# Setup analysis dict
analyse = Analysis.from_files(
    structure_path=os.path.join(directory, "CONTCAR"),
    icoxxlist_path=os.path.join(directory, "ICOHPLIST.lobster"),
    coxxcar_path=os.path.join(directory, "COHPCAR.lobster"),
    charge_path=os.path.join(directory, "CHARGE.lobster"),
    which_bonds="all",
    orbital_resolved=True,
    orbital_cutoff=0.05,
)

# Setup Description dict
describe = Description(analysis_object=analyse)
describe.write_description()



# Automatic interactive plots
describe.plot_interactive_cohps(ylim=[-10, 2], xlim=[-4, 4], orbital_resolved=True, label_resolved=False)

# different dicts that summarize the results
print(analyse.condensed_bonding_analysis)
print(analyse.final_dict_bonds)
print(analyse.final_dict_ions)
