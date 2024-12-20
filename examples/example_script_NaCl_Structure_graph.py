from lobsterpy.structuregraph.graph import LobsterGraph

graph_NaCl_all = LobsterGraph(
    path_to_poscar="./NaCl_comp_range/CONTCAR.gz",
    path_to_charge="./NaCl_comp_range/CHARGE.lobster.gz",
    path_to_cohpcar="./NaCl_comp_range/COHPCAR.lobster.gz",
    path_to_icohplist="./NaCl_comp_range/ICOHPLIST.lobster.gz",
    add_additional_data_sg=True,
    path_to_icooplist="./NaCl_comp_range/ICOOPLIST.lobster.gz",
    path_to_icobilist="./NaCl_comp_range/ICOBILIST.lobster.gz",
    path_to_madelung="./NaCl_comp_range/MadelungEnergies.lobster.gz",
    which_bonds="all",
    start=None,
)
print(graph_NaCl_all.sg)
