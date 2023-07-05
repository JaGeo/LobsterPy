import numpy as np
import matplotlib.pyplot as plt
import fnmatch
from lobsterpy.misc.constants import ATOM_TYPES


class DOS:
    """
    DOS class to hold all contents of the DOSCAR.lobster file as a dict. For each state (total DOS, orbital-wise DOSs),
    the dict contains all information, so there is some redudancy in the data TODO: remove redundancy
    Allows convenient plotting of DOS via simple string input.
    """

    def __init__(self, doscar_filename='DOSCAR.lobster'):
        self.dos_states = {}
        self._read_doscar(doscar_filename)

    def get_dos_states(self):
        return self.dos_states.keys()

    def get_dos_state(self, state):
        return self.dos_states[state]

    def print_plot(self, states_string, e_range=[]):
        max_dos_val = 0.0
        fig, ax = plt.subplots()
        fig.set_size_inches(8.0,12.0)
        for dos_state in states_string:
            this_dos = np.array([0.0 for i in range(self.dos_states['total_dos']['n_points'])])
            energy = self.dos_states['total_dos']['energy']
            for sub_state in dos_state.split("+"):
                matched_states = fnmatch.filter(self.dos_states.keys(), sub_state)
                if len(matched_states) == 0:
                    print("WARNING: No states matched for {state}".format(state=sub_state))
                for matched_state in matched_states:
                    this_dos += self.dos_states[matched_state]['dos']
            max_dos_val = max(max_dos_val, this_dos.max())
            ax.plot(this_dos, energy, label=dos_state)
        ax.set_xlabel('DOS')
        ax.set_ylabel("Energy (eV)")
        if e_range:
            ax.set_ylim(e_range)
        ax.set_xlim([0, max_dos_val/0.9])
        ax.axhline(y=0.0, color='black', linestyle='--')
        ax.axvline(x=0.0, color='black', linestyle='-')
        ax.text(ax.get_xlim()[1]/0.90, 0, "$\epsilon_\mathrm{F}$", color="black", ha="right", va="center")
        ax.legend(fontsize=14)
        fig.savefig("DOS.eps", format='eps', dpi=100, bbox_inches='tight')

    def _read_doscar(self, filename):
        """
        Reads the DOSCAR.lobster file and stores the information within a member variable "dos_states"
        """
        with open(filename, 'r') as doscar:
            # Skip header
            for i in range(5):
                doscar.readline()
            # Iterate through atom-wise DOSs
            atom_iterator = 1
            while True:
                this_line = doscar.readline()
                # Stop loop when we reached the end of the file
                if not this_line:
                    break
                # Get information about this atom from the first line of the block
                info_segments = this_line.split(';')
                if len(info_segments) == 1:
                    # If there is no information about atom, then this is the total DOS
                    this_state = ['total_dos']
                else:
                    # Else, we can extract atomic number, infer atomic type and extract information about atomic orbitals
                    atomic_number = int(info_segments[1].split()[1])
                    atom_type = ATOM_TYPES[atomic_number]
                    this_state = ["{atom_type}.{atom_iterator}.{atomic_orbital}".format(atom_type=atom_type,
                                                                                               atom_iterator=atom_iterator,
                                                                                               atomic_orbital=atomic_orbital)
                                         for atomic_orbital in info_segments[2].split()]
                    atom_iterator += 1
                energy_infos = info_segments[0].split()
                e_min = float(energy_infos[0])
                e_max = float(energy_infos[1])
                e_f = float(energy_infos[3])
                n_points = int(energy_infos[2])
                # Prepare arrays for energy range and all atomic orbital DOSs on atom
                energy = []
                state_dos = [[] for i in range(len(this_state))]
                # Loop through list of DOS points and store in arrays
                for i in range(n_points):
                    dos_line = [float(val) for val in doscar.readline().split()]
                    energy.append(dos_line[0])
                    for state in range(len(state_dos)):
                        state_dos[state].append(dos_line[state + 1])
                # Store in dictionary
                for i, state in enumerate(this_state):
                    self.dos_states[state] = {
                        'energy': np.array(energy),
                        'dos': np.array(state_dos[i]),
                        'e_min': e_min,
                        'e_max': e_max,
                        'e_f': e_f,
                        'n_points': n_points,
                        'atom_type': state.split('.')[0],
                        'shell': state.split('.')[-1][:2],
                        'atomic_orbital': state.split('.')[-1]
                    }
