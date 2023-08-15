import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.api.types import CategoricalDtype
from tqdm import tqdm

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks")
sns.set(font="Helvetica")

plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams.update({"font.size": 7})


chr_order = CategoricalDtype(
    [
        "CHR_1",
        "CHR_2",
        "CHR_3",
        "CHR_4",
        "CHR_5",
        "CHR_6",
        "CHR_7",
        "CHR_8",
        "CHR_9",
        "CHR_10",
        "CHR_11",
        "CHR_12",
        "CHR_13",
        "CHR_14",
        "CHR_15",
        "CHR_16",
        "CHR_MT",
    ],
    ordered=True,
)


def get_chr_limits(snp_position):
    #### Delimitations of chromosomes based on SNPs index
    chr_limits = []
    for chr in np.unique(snp_position["chrom"].values):
        TMP = snp_position[snp_position["chrom"] == chr]
        chr_limits.append(
            [
                chr,
                int(TMP["SNP"].describe()["min"]),
                int(TMP["SNP"].describe()["max"]),
                int(TMP["position"].describe()["min"]),
                int(TMP["position"].describe()["max"]),
            ]
        )
    chr_limit_table = pd.DataFrame(
        chr_limits,
        columns=[
            "CHR",
            "SNP_IDX_FIRST",
            "SNP_IDX_LAST",
            "SNP_POS_FIRST",
            "SNP_POS_LAST",
        ],
    ).sort_values("SNP_IDX_FIRST")
    chr_limit_table["CHR"] = chr_limit_table["CHR"].astype(chr_order)
    print(chr_limit_table)
    return chr_limit_table


def create_genotype_map(genotype_matrix):
    #### Remove snp_id column
    cols = genotype_matrix.columns
    genotype_matrix = genotype_matrix[cols[1:]]
    print(genotype_matrix.shape)
    f, ax = plt.subplots(nrows=1, ncols=1)
    ax.grid(False)
    plt.imshow(
        genotype_matrix.T.values,
        vmin=-1,
        vmax=1,
        extent=[1, 12054, 357, 1],
        cmap="viridis",
        interpolation="None",
        aspect="auto",
    )
    plt.colorbar()

    # ToDo: modify output paths
    f.savefig("../manuscript/figures/EXT_FIGURE_1/FIGURE_1A.png", dpi=300)
    f.savefig("../manuscript/figures/EXT_FIGURE_1/FIGURE_1A.eps", dpi=300)


def format_genotype(x):
    if x == -1:
        return 0
    elif x == 1:
        return 2
    elif x == 0:
        return -1


def linkage_disequilibrium(genotype_matrix):
    genotype_matrix = genotype_matrix.replace(0, "NA")
    print(genotype_matrix)

    f = open(
        "../data/genotype_information/piQTL_genotype_matrix_for_LD_dec2022.txt", "w"
    )
    for idx in tqdm(genotype_matrix.index):
        for col in genotype_matrix.columns[1:]:
            line = "\t".join(
                [
                    str(genotype_matrix["snp_id"][idx]),
                    col,
                    str(genotype_matrix.loc[idx, col]),
                    str(1),
                ]
            )
            f.write(line + "\n")
    f.close()


if __name__ == "__main__":
    genotype_matrix = pd.read_csv(
        "../data/genotype_information/piQTL_genotype_matrix_dec2022.txt"
    )
    snp_position = pd.read_csv(
        "../data/genotype_information/snps_annotations_genome-version-3-64-1.txt"
    ).rename(columns={"snp_id": "SNP"})
    print(genotype_matrix.shape)
    chr_limits = get_chr_limits(snp_position)
    print("Save panel A for Figure 1")
    create_genotype_map(genotype_matrix)
    print("Input for generating panel B for Figure 1 with Rscript calculate_LD.R ")
