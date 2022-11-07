import os
import random
import warnings

import numpy as np
import plotly.graph_objs as go

from lobsterpy.structuregraph.graph import LobsterGraph
from pymatgen.electronic_structure.cohp import CompleteCohp
from itertools import product, permutations

warnings.filterwarnings(action='ignore')



def get_coord_transformation_matrix(
    a: float, b: float, c: float,
    alpha: float, beta: float, gamma: float, is_rad: bool = False
) -> np.array:
    """
    function to get matrix which converts fractional to cartesian coordinates in a given crystal system

    :param a:
    :param b:
    :param c:
    :param alpha:
    :param beta:
    :param gamma:
    :param is_rad: True if the angles are in radians, False if in degree (default: False)
    :return: transformation matrix as numpy array
    """

    if not is_rad:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    a_vec = np.array([a, 0, 0])
    b_vec = np.array([0, b, 0])
    c_vec = np.array([0, 0, c])

    x = np.array([a, 0, 0])
    y = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    z = np.array([c * np.cos(beta), c * np.cos(alpha), c * np.sin(gamma)])
    # https://en.wikipedia.org/wiki/Fractional_coordinates
    # spat = np.dot(a_vec, np.cross(b_vec, c_vec))
    #z = np.array(
    #    [
    #        c * np.cos(beta),
    #        c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)),
    #        spat / (a * b * np.sin(gamma))
    #    ]
    #)

    return np.stack((x, y, z), axis=-1)



def get_cell_axes(
    a: float, b: float, c: float,
    alpha: float, beta: float, gamma: float, is_rad: bool = False
) -> list:
    axes = []

    if not is_rad:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    # https://en.wikipedia.org/wiki/Fractional_coordinates
    #a_vec = np.array([a, 0, 0])
    #b_vec = np.array([0, b, 0])
    #c_vec = np.array([0, 0, c])
    #spat = np.dot(a_vec, np.cross(b_vec, c_vec))
    #x = np.array([a, 0, 0])
    #y = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    #z = np.array(
    #    [
    #        c * np.cos(beta),
    #        c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)),
    #        spat / (a * b * np.sin(gamma))
    #    ]
    #)

    origin = np.array([0, 0, 0])
    # cartesian x-axis
    x = np.array([a, 0, 0])
    axes.append((origin, x))
    # cartesian y-axis
    y = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
    axes.append((origin, y))
    # cartesian z-axis
    z = np.array([c * np.cos(beta), c * np.cos(alpha), c * np.sin(gamma)])
    axes.append((origin, z))
    # cartesian x-axis parallel in y-direction
    axes.append((y, y + x))
    # cartesian x-axis parallel in z-direction
    axes.append((z, z + x))
    # cartesian x-axis parallel in yz-direction
    axes.append((y + z, y + z + x))
    # cartesian y-axis parallel in x-direction
    axes.append((x, x + y))
    # cartesian y-axis parallel in z-direction
    axes.append((z, z + y))
    # cartesian y-axis parallel in xz-direction
    axes.append((x + z, x + z + y))
    # cartesian z-axis parallel in x-direction
    axes.append((x, x + z))
    # cartesian z-axis parallel in y-direction
    axes.append((y, y + z))
    # cartesian z-axis parallel in xy-direction
    axes.append((x + y, x + y + z))

    return axes



def get_primitive_cell(lobstergraph: LobsterGraph, completecohp: CompleteCohp) -> (list, dict, dict):
    """
    function for building first full primitive cell with edges

    :param lobstergraph: graph object containing information about nodes (atoms) and edges (bonds)
    :return cells: list of primitive cells, contains one cell after this function is finished
    :return atoms: dictionary of atoms in primitive cell, keys are enumerated atom numbers, values element symbol and
                   element number
    :return eq_atoms: dictionary of atoms equivalent the already existing ones in the graph object, keys are
                      enumerated atoms numbers, values are enumerated atom numbers of equivalent new atoms
    """

    # initialization block
    cell = {
        "atoms": {},
        "axes": [],
        "edges": {},
    }
    structure = lobstergraph.sg.structure
    num_atoms = len(structure.frac_coords)
    eq_atoms = dict()
    tol = 0.01

    # iterate over all already existing nodes, save cartesian coordinates, fractional coordinates, element and
    # element number
    for i, frac_coord in enumerate(structure.frac_coords.copy()):
        cell["atoms"][i] = {
            "element": str(structure[i].specie),
            "number": structure[i].specie.number,
            "frac_coord": frac_coord,
        }
        eq_atoms[i] = []
        # for every node find fractional coordinates with entries 0 or 1 in some dimension, for these nodes equivalent
        # nodes/coordinates (due to translational symmetry) will be created
        indices0 = []
        indices1 = []
        # check the comparison
        for j, c in enumerate(frac_coord):
            if c < 0.01:
                indices0.append(j)
            elif c > 0.99:
                indices1.append(j)
        n0 = len(indices0)
        n1 = len(indices1)
        add0 = [np.zeros(shape=3)]
        add1 = [np.zeros(shape=3)]
        # if there is more than one entry with 0 or 1 create all possible permutations for the translational shift
        for j in range(1, n0 + 1):
            v = [1] * j + [0] * (n0 - n1 - j)
            ps0 = set(permutations(v))
            for p0 in ps0:
                a0 = np.zeros(shape=3)
                for index, permutation in zip(indices0, p0):
                    a0[index] = permutation
                add0.append(a0)
        for j in range(1, n1 + 1):
            v = [-1] * j + [0] * (n1 - n0 - j)
            ps1 = set(permutations(v))
            for p1 in ps1:
                a1 = np.zeros(shape=3)
                for index, permutation in zip(indices1, p1):
                    a1[index] = permutation
                add1.append(a1)
        adds = product(add0, add1)
        # create new nodes/coordinates by adding translational shift vectors
        # save new fractional coordinates, cartesian coordinates, element and element number
        for a0, a1 in adds:
            if np.linalg.norm(a0 + a1) > 0:
                cell["atoms"][num_atoms] = {
                    "element": str(structure[i].specie),
                    "number": structure[i].specie.number,
                    "frac_coord": frac_coord + a0 + a1,
                }
                eq_atoms[i].append(num_atoms)
                num_atoms += 1

    # iterate over all already existing nodes
    for node1, node2, data in lobstergraph.sg.graph.edges.data():
        # extract edge properties
        bond_label = data["bond_label"]
        cohp_data = completecohp.get_cohp_by_label(label=bond_label).as_dict()
        spinup_cohps = cohp_data["COHP"]["1"]
        spindown_cohps = cohp_data["COHP"]["-1"]
        energies = cohp_data["energies"]
        fermi_energy = cohp_data["efermi"]
        x = [-(spinup_cohps[j] + spindown_cohps[j]) for j, _ in enumerate(spinup_cohps)]
        y = [energies[j] - fermi_energy for j, _ in enumerate(energies)]
        frac_coord1 = cell["atoms"][node1]["frac_coord"]
        frac_coord2 = cell["atoms"][node2]["frac_coord"] + data["to_jimage"]
        # check if edges lies within cell
        if (-tol <= frac_coord2[0] <= 1+tol) and \
           (-tol <= frac_coord2[1] <= 1+tol) and \
           (-tol <= frac_coord2[2] <= 1+tol):
            cell["edges"][bond_label] = {
                "frac_coords": [(frac_coord1, frac_coord2)],
                "cohp_plot": (x, y),
                "bond_length": data["bond_length"],
                "icobi": data["ICOBI"],
                "icoop": data["ICOOP"],
                "icohp": data["ICOHP"],
                "icohp_bonding_perc": data["ICOHP_bonding_perc"]
            }
        # iterate over all edges that are equivalent to the current one due to translational symmetry
        for eq_atom in eq_atoms[node1]:
            start = cell["atoms"][eq_atom]["frac_coord"]
            shift = start - frac_coord1
            end = frac_coord2 + shift
            # check if edges lies within cell
            if (-tol <= end[0] <= 1+tol) and \
               (-tol <= end[1] <= 1+tol) and \
               (-tol <= end[2] <= 1+tol):
                # if statement more elegant !
                try:
                    cell["edges"][bond_label]["frac_coords"].append((start, end))
                except:
                    cell["edges"][bond_label] = {
                        "frac_coords": [(start, end)],
                        "cohp_plot": (x, y),
                        "bond_length": data["bond_length"],
                        "icobi": data["ICOBI"],
                        "icoop": data["ICOOP"],
                        "icohp": data["ICOHP"],
                        "icohp_bonding_perc": data["ICOHP_bonding_perc"]
                    }

    return cell



def get_primitive_supercell(
        lobstergraph: LobsterGraph, cell: dict, frac_to_cart_matrix: np.ndarray
) -> dict:
    """
    function to build a primitive supercell from a primitive cell based on connectivity information in graph
    object ("to_jimage" vectors)

    :param lobstergraph: graph object containing information about nodes (atoms) and edges (bonds)
    :param cells: list containing primitive cell to duplicate in order to build supercell
    :param atoms: dictionary of atoms in primitive cell, keys are enumerated atom numbers, values element symbol and
                  element number
    :param frac_to_cart_matrix: matrix to convert fractional coordinates to cartesian coordinates
    :return: cells: list containing all primitive cells that make up primitive supercell
    """

    data = list(lobstergraph.sg.graph.edges.data())

    xs = [vec["to_jimage"][0] for _, _, vec in data]
    # maximum x-component of all "to_jimage" vectors determines size of supercell in +x direction
    x_min = min(xs)
    # minimum x-component of all "to_jimage" vectors determines size of supercell in -x direction
    x_max = max(xs)

    ys = [vec["to_jimage"][1] for _, _, vec in data]
    # maximum y-component of all "to_jimage" vectors determines size of supercell in +y direction
    y_min = min(ys)
    # minimum y-component of all "to_jimage" vectors determines size of supercell in -y direction
    y_max = max(ys)

    zs = [vec["to_jimage"][2] for _, _, vec in data]
    # maximum z-component of all "to_jimage" vectors determines size of supercell in +z direction
    z_min = min(zs)
    # minimum z-component of all "to_jimage" vectors determines size of supercell in -z direction
    z_max = max(zs)

    # get cell parameters of supercell (necessary for building extended supercell)
    cell["a"] = 1 + abs(x_max) + abs(x_min)
    cell["b"] = 1 + abs(y_max) + abs(y_min)
    cell["c"] = 1 + abs(z_max) + abs(z_min)

    # iterate over x, y, z dimension
    num_atoms = len(cell["atoms"])
    for i, (dim_min, dim_max) in enumerate([(x_min, x_max), (y_min, y_max), (z_min, z_max)]):
        # create new cell by shifting existing one in x, y or z direction
        shift = np.array([0, 0, 0])
        shift[i] = 1.0
        new_atoms = {}
        new_axes = []
        new_edges = {}
        # repeat every cell x_min/y_min/z_min times in negative direction and x_max/y_max/z_max times in
        # positive direction
        for j in [k for k in range(dim_min, dim_max + 1) if k != 0]:
            # calculate coordinates of new nodes (atoms)
            for atom in cell["atoms"].values():
                new_atoms[num_atoms] = {
                    "element": atom["element"],
                    "number": atom["number"],
                    "frac_coord": atom["frac_coord"] + j * shift,
                }
                num_atoms += 1
            # calculate coordindates of new axes
            for start, end in cell["axes"]:
                new_start = start + np.dot(frac_to_cart_matrix, j * shift)
                new_end = end + np.dot(frac_to_cart_matrix, j * shift)
                new_axes.append((new_start, new_end))
            # calculate coordinates of new edges (bonds), add them to bond labels based on equivalence
            for bond_label, edge in cell["edges"].items():
                new_edges[bond_label] = []
                for start, end in edge["frac_coords"]:
                    new_start = start + j * shift
                    new_end = end + j * shift
                    new_edges[bond_label].append((new_start, new_end))
            # add new nodes to cell
            cell["atoms"] = {**cell["atoms"], **new_atoms}
            # add new axes to cell
            cell["axes"] += new_axes
            # add new edges by bond label to cell
            for bond_label in new_edges.keys():
                for frac_coords in new_edges[bond_label]:
                    cell["edges"][bond_label]["frac_coords"].append(frac_coords)

    return cell



def get_extended_primitive_supercell(cell: dict, vecs: list, frac_to_cart_matrix: np.ndarray) -> dict:
    """
    function to extend the computed primitive supercell in desired spatial directions

    :param cell: dictionary defining the primitive supercell to be extended
    :param vecs: list of shift vectors that define how often and in which direction the primitive supercell is
                 repeated, format [n, a, b, c] where n is how often the cell should be repeated and
                 a, b, c = -1, 0, 1 define in which direction the cell will be repeated
    :param frac_to_cart_matrix: matrix to convert fractional coordinates to cartesian coordinates in the
                                given crystal system
    :return: cell: dictionary defining the extended primitive supercell
    """

    # initialization block
    num_atoms = len(cell["atoms"])
    new_atoms = {}
    new_axes = []
    new_edges = {bond_label: [] for bond_label, data in cell["edges"].items()}

    # apply all user defined shift vectors
    for vec in vecs:
        # number of times to repeat cell in certain direction
        n = vec[0]
        # direction in which the cell will be repeated
        shift = np.array(vec[1:])
        shift[0] *= cell["a"]
        shift[1] *= cell["b"]
        shift[2] *= cell["c"]
        # repeat the cell 1, 2, ..., n times in the certain direction
        for i in range(1, n+1):
            # calculate coordinates of new nodes (atoms)
            for atom in cell["atoms"].values():
                new_atoms[num_atoms] = {
                    "element": atom["element"],
                    "number": atom["number"],
                    "frac_coord": atom["frac_coord"] + i * shift
                }
                num_atoms += 1
            # calculate coordindates of new axes
            for start, end in cell["axes"]:
                new_start = start + np.dot(frac_to_cart_matrix, i * shift)
                new_end = end + np.dot(frac_to_cart_matrix, i * shift)
                new_axes.append((new_start, new_end))
            # calculate coordinates of new edges (bonds), add them to bond labels based on equivalence
            for bond_label, edge in cell["edges"].items():
                for start, end in edge["frac_coords"]:
                    new_start = start + i * shift
                    new_end = end + i * shift
                    new_edges[bond_label].append((new_start, new_end))

    # add new nodes to cell
    cell["atoms"] = {**cell["atoms"], **new_atoms}
    # add new axes to cell
    cell["axes"] += new_axes
    # add new edges by bond label to cell
    for bond_label in new_edges.keys():
        for frac_coords in new_edges[bond_label]:
            cell["edges"][bond_label]["frac_coords"].append(frac_coords)

    return cell



def get_structure_plot(
        path_to_poscar: str,
        path_to_charge: str,
        path_to_icobilist: str,
        path_to_icooplist: str,
        path_to_icohplist: str,
        path_to_cohpcar: str,
        path_to_madelung: str,
        vecs: list = None
) -> go.Figure:
    """
    Creation of an interactive 3D plot of a compound's primitive supercell, containing information about site and
    bond properties.

    :param lobstergraph: LobsterGraph object, contains information about connectivity, site and bond properties in a
                         graph-like manner
    :return: fig: visualization of primitive supercell by 3D scatter plot
    """

    # initialization block
    lobstergraph = LobsterGraph(
        path_to_poscar=path_to_poscar,
        path_to_charge=path_to_charge,
        path_to_icobilist=path_to_icobilist,
        path_to_icooplist=path_to_icooplist,
        path_to_icohplist=path_to_icohplist,
        path_to_cohpcar=path_to_cohpcar,
        path_to_madelung=path_to_madelung,
        which_bonds="all",
        # start=-2,
        add_additional_data_sg=True
    )

    # already included in lobstergraph !
    completecohp = CompleteCohp.from_file(
        fmt="LOBSTER", filename=path_to_cohpcar, structure_file=path_to_poscar
    )

    atom_number = []

    node_x = []
    node_y = []
    node_z = []

    axis_x = []
    axis_y = []
    axis_z = []

    # get matrix to transform fractional to cartesian coordinates
    a = lobstergraph.sg.structure.lattice.a
    b = lobstergraph.sg.structure.lattice.b
    c = lobstergraph.sg.structure.lattice.c
    alpha = lobstergraph.sg.structure.lattice.alpha
    beta = lobstergraph.sg.structure.lattice.beta
    gamma = lobstergraph.sg.structure.lattice.gamma
    # use function from pymatgen !
    frac_to_cart_matrix = get_coord_transformation_matrix(a, b, c, alpha, beta, gamma)

    # build primitive cell
    cell = get_primitive_cell(lobstergraph, completecohp)

    # use function from pymatgen !
    # add cell axes to primitive cell
    cell["axes"] = get_cell_axes(a, b, c, alpha, beta, gamma)

    # build primitive supercell from primitive cell
    cell = get_primitive_supercell(lobstergraph, cell, frac_to_cart_matrix)

    # build extended primitive supercell from primitive supercell
    if vecs:
        cell = get_extended_primitive_supercell(cell=cell, vecs=vecs, frac_to_cart_matrix=frac_to_cart_matrix)

    # layout of plot axes
    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
        showspikes=False
    )

    # overall layout of plot
    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(axis),
            yaxis=dict(axis),
            zaxis=dict(axis),
        ),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        hovermode="closest",
        height=820,
        width=1195
    )

    fig = go.Figure(layout=layout)

    # add nodes (atoms) with color depending on element to plot
    for atom in cell["atoms"].values():
        coord = np.dot(frac_to_cart_matrix, atom["frac_coord"])
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_z.append(coord[2])
        atom_number.append(atom["number"])
    node_trace = go.Scatter3d(
        x=node_x,
        y=node_y,
        z=node_z,
        mode="markers",
        hoverinfo="none",
        marker=dict(
            symbol="circle",
            size=6,
            color=atom_number,
            colorscale="Viridis",
            line=dict(color="rgb(50,50,50)", width=0.5),
        ),
    )
    fig.add_trace(node_trace)

    # add axes of primitive cells to plot
    for start, end in cell["axes"]:
        axis_x += [start[0], end[0], None]
        axis_y += [start[1], end[1], None]
        axis_z += [start[2], end[2], None]
    axes_trace = go.Scatter3d(
        x=axis_x,
        y=axis_y,
        z=axis_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="grey", width=1)
    )
    fig.add_trace(axes_trace)

    # add edges with edge properties as custom hover data as separate traces to plot, to make them separately
    # accessible for hover events
    for edge in cell["edges"].values():
        for start, end in edge["frac_coords"]:
            start = np.dot(frac_to_cart_matrix, start)
            end = np.dot(frac_to_cart_matrix, end)
            fig.add_trace(
                go.Scatter3d(
                    x=[start[0], end[0], None],
                    y=[start[1], end[1], None],
                    z=[start[2], end[2], None],
                    mode="lines",
                    line={
                        "width": 2,
                        "color": "black"
                    },
                    hoverinfo="none",
                    customdata=[
                        [
                            edge["cohp_plot"],
                            edge["bond_length"],
                            edge["icobi"],
                            edge["icoop"],
                            edge["icohp"],
                            edge["icohp_bonding_perc"]
                        ],
                        [
                            edge["cohp_plot"],
                            edge["bond_length"],
                            edge["icobi"],
                            edge["icoop"],
                            edge["icohp"],
                            edge["icohp_bonding_perc"]
                        ],
                        None
                    ]
                )
            )

    return fig


# TODO
def create_graph_plot(lobstergraph: LobsterGraph):
    atom_number = []

    node_x = []
    node_y = []
    node_z = []

    edge_x = []
    edge_y = []
    edge_z = []

    frac_coords = []
    for i, _ in enumerate(lobstergraph.sg.structure.frac_coords):
        new_frac_coord = [(-1)**i * i, 0, 0]
        frac_coords.append(new_frac_coord)
        atom_number.append(lobstergraph.sg.structure[i].specie.number)

    edges = []
    for i, (node1, node2, _) in enumerate(lobstergraph.sg.graph.edges.data()):
        start = frac_coords[node1][0]
        stop = frac_coords[node2][0]
        mid = 0.5 * (stop - start)
        x = np.linspace(start, stop)
        y = (x - (start+mid))**2 - mid**2
        z = [0] * len(x)
        edges.append((x, y, z))

    # layout of plot axes
    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
        showspikes=False
    )

    # overall layout of plot
    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(axis),
            yaxis=dict(axis),
            zaxis=dict(axis),
        ),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        hovermode="closest",
        height=820,
        width=1195
    )

    fig = go.Figure(layout=layout)

    for frac_coord, number in zip(frac_coords, atom_number):
        node_x.append(frac_coord[0])
        node_y.append(frac_coord[1])
        node_z.append(frac_coord[2])
        atom_number.append(number)
    node_trace = go.Scatter3d(
        x=node_x,
        y=node_y,
        z=node_z,
        mode="markers",
        hoverinfo="none",
        marker=dict(
            symbol="circle",
            size=6,
            color=atom_number,
            colorscale="Viridis",
            line=dict(color="rgb(50,50,50)", width=0.5),
        ),
    )
    fig.add_trace(node_trace)

    #edge_trace = go.Scatter3d(
    #    x=[-5, 5, None],
    #    y=[0, 0, None],
    #    z=[0, 0, None],
    #    mode="lines",
    #    customdata=None,
    #    hoverinfo="none",
    #    line=dict(color="black", width=2),
    #)

    for edge in edges:
        n = len(edge)
        for i in range(1, n):
            start = edge[i-1]
            end = edge[i]
            edge_x += [start[0], end[0], None]
            edge_y += [start[1], end[1], None]
            edge_z += [start[2], end[2], None]
    edge_trace = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode="lines",
        customdata=None,
        hoverinfo="none",
        line=dict(color="black", width=2),
    )
    fig.add_trace(edge_trace)

    #fig.add_trace(
    #    go.Scatter3d(

    #    )
    #)

    fig.show()



def get_dummy_cohp_plot() -> go.Figure:
    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        visible=False,
        title="",
        showspikes=False
    )

    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=axis,
            yaxis=axis,
        ),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        height=440,
        width=600,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False)
    )

    return go.Figure(layout=layout)



def get_random_structure() -> str:
    filepath = os.path.expanduser("~/automationresults")
    mps = os.listdir(filepath)
    rand_index = random.randrange(len(mps))
    dir = mps[rand_index]
    path = os.path.join(filepath, dir)
    print(f"RANDOMLY CHOSEN STRUCTURE: {dir}")
    return path



def get_chosen_structure(file: str = None) -> str:
    if file is None:
        path = get_random_structure()
    else:
        path = os.path.join(os.path.expanduser("~/automationresults"), file)
    return path
