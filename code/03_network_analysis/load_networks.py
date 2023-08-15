import numpy as np
import pandas as pd

from tqdm import tqdm


def format_col(x):
    x = x.split(".")[1]
    return x


def get_string_network(root_path):
    string_network = pd.read_csv(
        f"{root_path}networks/STRING/4932.protein.physical.links.detailed.v11.5.txt",  # ToDo: no input file
        sep=" ",
    )
    string_network["protein1"] = string_network["protein1"].apply(format_col)
    string_network["protein2"] = string_network["protein2"].apply(format_col)
    string_network["experimental"] = string_network["experimental"] / 1000

    final_string_network = string_network[
        string_network["experimental"] > 0.5
    ].reset_index(drop=True)
    return final_string_network[["protein1", "protein2", "experimental"]]


def get_biogrid_network(root_path, exp_type):
    biogrid_network = pd.read_table(
        f"{root_path}networks/BioGrid/BIOGRID-ORGANISM-4.4.200.tab3/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.200.tab3.txt",  # ToDo: no input file
        low_memory=False,
    )
    biogrid_network = biogrid_network[
        (
            biogrid_network["Organism Name Interactor A"]
            == "Saccharomyces cerevisiae (S288c)"
        )
        & (
            biogrid_network["Organism Name Interactor B"]
            == "Saccharomyces cerevisiae (S288c)"
        )
    ]
    biogrid_network_exp = (
        biogrid_network[biogrid_network["Experimental System"] == exp_type]
        .reset_index(drop=True)
        .rename(
            columns={
                "Systematic Name Interactor A": "protein1",
                "Systematic Name Interactor B": "protein2",
            }
        )
    )
    return biogrid_network_exp[["protein1", "protein2"]]


def get_tarassov_intensity(intensity_matrix):
    #### Selecting edges of the network based on a given threshold
    edges = []
    for row in tqdm(intensity_matrix.index):
        for col in intensity_matrix.columns:
            intensity_val = intensity_matrix.at[row, col]
            edges.append([row, col, intensity_val])
    return edges


def get_tarassov_network(root_path, protein_abundance):
    index_ids = pd.read_csv(
        f"{root_path}networks/tarassov_updated_dataset/alpha_array.csv",  # ToDo: no input file
        names=["index"],
    )["index"]
    columns_ids = pd.read_csv(
        f"{root_path}networks/tarassov_updated_dataset/allInteractionsNames.csv",  # ToDo: no input file
        names=["columns"],
    )["columns"]
    intensity_matrix = pd.read_table(
        f"{root_path}networks/tarassov_updated_dataset/allInteractionsIntensities.txt",  # ToDo: no input file
        header=None,
    )

    # print(len(index_ids), len(columns_ids))
    # print(intensity_matrix.shape)

    intensity_matrix.index = index_ids
    intensity_matrix.columns = columns_ids

    true_idx_pa = []
    for idx in tqdm(intensity_matrix.index):
        if idx in protein_abundance["Systematic Name"].values:
            if idx not in true_idx_pa:
                true_idx_pa.append(idx)
    true_col_pa = []
    for col in tqdm(intensity_matrix.columns):
        if col in protein_abundance["Systematic Name"].values:
            if col not in true_col_pa:
                true_col_pa.append(col)

    ### Restrict table to protein with available abundance
    intensity_matrix = intensity_matrix[intensity_matrix.index.isin(true_idx_pa)][
        true_col_pa
    ]
    print(intensity_matrix.shape)

    tarassov_edges = get_tarassov_intensity(intensity_matrix)

    tarassov_network = pd.DataFrame(
        tarassov_edges, columns=["protein1", "protein2", "intensity"]
    )
    return tarassov_network


def get_cellmap_edge_list(similarity_matrix, threshold):
    #### Selecting edges of the network based on a given threshold
    edges = []
    for i in tqdm(range(len(similarity_matrix.values))):
        row = similarity_matrix.values[i]
        nA = similarity_matrix.columns[i]
        for j in range(len(row)):
            indice = row[j]
            nB = similarity_matrix.columns[j]
            if indice > threshold:
                if nA != nB:
                    edges.append((nA, nB, indice))
    return edges


def get_cellmap_network(root_path):
    cellMap_similarity_matrix = pd.read_table(
        f"{root_path}/networks/cellMap/similarity_matrices/cc_ALL.txt",  # ToDo: no input file
        sep="\t",
        low_memory=False,
    )
    cellMap_similarity_matrix = cellMap_similarity_matrix.drop(
        columns=["Unnamed: 0", "Unnamed: 1"]
    )
    cellMap_similarity_matrix.columns = cellMap_similarity_matrix.values[0]
    cellMap_similarity_matrix = cellMap_similarity_matrix.drop(0).reset_index(drop=True)

    cellMap_similarity_matrix = cellMap_similarity_matrix.astype(float)
    cellMap_similarity_matrix = cellMap_similarity_matrix.replace(np.nan, -1)

    edges = get_cellmap_edge_list(cellMap_similarity_matrix, 0.2)

    cellMap_network = pd.DataFrame(edges, columns=["protein1", "protein2", "pcc"])
    return cellMap_network
