import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_scatter_density
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from pandas.api.types import CategoricalDtype
from scipy.stats import gaussian_kde

# set matplotlib default parameters
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = "Helvetica"
mpl.rcParams["font.size"] = 5
mpl.rcParams["axes.titlesize"] = 6
mpl.rcParams["xtick.labelsize"] = 5
mpl.rcParams["xtick.major.size"] = 2
mpl.rcParams["xtick.minor.size"] = 1.2
mpl.rcParams["ytick.labelsize"] = 5
mpl.rcParams["ytick.major.size"] = 2
mpl.rcParams["ytick.minor.size"] = 1.2
mpl.rcParams["lines.linewidth"] = 1.2
mpl.rcParams["lines.markersize"] = 2
mpl.rcParams["lines.markeredgewidth"] = 0.5
mpl.rcParams["boxplot.flierprops.markersize"] = 2
mpl.rcParams["boxplot.flierprops.markeredgewidth"] = 0.5
mpl.rcParams["legend.fontsize"] = 5

CM = 1 / 2.54  # centimeters in inches
SINGLE = 12 * CM
DOUBLE = 18 * CM

DRUG_ORDER = CategoricalDtype(
    ["noDrug", "5-FC", "Fluconazole", "Metformin", "Trifluoperazine"], ordered=True
)

sns.set_context("paper", font_scale=1)


def fitness_comparison(PPI, DRUG, DOUBLE):
    f, axes = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(DOUBLE, 12 * CM))
    i = 0
    for MTX in ["noMTX", "MTX"]:
        REPS = []
        LT = []
        for REP in ["A", "B"]:
            PPI_table = pd.read_csv(
                f"../../results/01_barcode_count/per_PPI/{PPI}_read_number.csv"
            )
            R = PPI_table[
                ["strain_id", f"T0_noMTX_noDrug.{PPI}", f"{REP}_{MTX}_{DRUG}.{PPI}"]
            ]
            R[f"{REP}_{MTX}_{DRUG}.{PPI}"] = R[f"{REP}_{MTX}_{DRUG}.{PPI}"] + 1
            R[f"T0.{PPI}"] = R[f"T0_noMTX_noDrug.{PPI}"] + 1

            R[f"{REP}_{MTX}_{DRUG}.{PPI}_RPM"] = (
                R[f"{REP}_{MTX}_{DRUG}.{PPI}"] / sum(R[f"{REP}_{MTX}_{DRUG}.{PPI}"])
            ) * 10**6
            R[f"T0.{PPI}_RPM"] = (R[f"T0.{PPI}"] / sum(R[f"T0.{PPI}"])) * 10**6
            R["logratio_Fitness"] = np.log2(
                R[f"{REP}_{MTX}_{DRUG}.{PPI}_RPM"] / R[f"T0.{PPI}_RPM"]
            )
            RPM_cols = [col for col in R.columns if "RPM" in col]
            RPM_cols.insert(0, "strain_id")
            RPM_cols.insert(-1, "logratio_Fitness")
            REPS.append(R[RPM_cols])

        #### Fitness comparison
        REP_MERGED = pd.merge(
            REPS[0], REPS[1], on=["strain_id", f"T0.{PPI}_RPM"], suffixes=["_A", "_B"]
        )
        REP_MERGED = REP_MERGED[
            REP_MERGED["strain_id"].isin(
                [
                    "17",
                    "17_2",
                    "40",
                    "42_2",
                    "180",
                    "180_2",
                    "599_ref1",
                    "599_ref2",
                    "43_ref1",
                    "43_ref2",
                ]
            )
        ]
        sns.scatterplot(
            data=REP_MERGED,
            x="logratio_Fitness_A",
            y="logratio_Fitness_B",
            linewidth=0,
            color="black",
            ax=axes[i],
        )
        axes[i].set_xlabel("$Fitness_{~Replicate~A}$")
        axes[i].set_ylabel("$Fitness_{~Replicate~B}$")
        sns.scatterplot(
            data=REP_MERGED[REP_MERGED["strain_id"].str.contains("_ref")],
            x="logratio_Fitness_A",
            y="logratio_Fitness_B",
            linewidth=0,
            color="red",
            ax=axes[i],
        )
        axes[i].axline([-1, -1], [1, 1], linestyle="dashed", color="grey")
        corr, pval = stats.pearsonr(
            REP_MERGED["logratio_Fitness_A"], REP_MERGED["logratio_Fitness_B"]
        )
        axes[i].annotate(
            f"Pearson's corr: {round(corr,3)}\npval: {'{:.2e}'.format(pval)}",
            xy=[0, -55],
            xycoords="axes points",
        )
        axes[i].set_title(f"{PPI} under {DRUG} ({MTX})".replace("_", ":"))
        axes[i].set_box_aspect(1)
        axes[i].yaxis.set_tick_params(which="both", labelleft=True)
        i += 1

    plt.show()


def using_mpl_scatter_density(fig, x, y):
    white_viridis = LinearSegmentedColormap.from_list(
        "white_viridis",
        [
            (0, "#ffffff"),
            (1e-20, "#440053"),
            (0.2, "#404388"),
            (0.4, "#2a788e"),
            (0.6, "#21a784"),
            (0.8, "#78d151"),
            (1, "#fde624"),
        ],
        N=256,
    )

    ax = fig.add_subplot(1, 1, 1, projection="scatter_density")
    density = ax.scatter_density(x, y, cmap=white_viridis)
    fig.colorbar(density, label="Number of points per pixel")

    return fig, ax


def display_density_plot(all_A_vs_B, MTX, SINGLE, root):
    f = plt.figure(figsize=(SINGLE, 9 * CM))

    all_A_vs_B = all_A_vs_B[~all_A_vs_B["strain_id"].str.contains("_ref")]
    xy = np.vstack([all_A_vs_B["logratio_Fitness_A"], all_A_vs_B["logratio_Fitness_B"]])
    z = gaussian_kde(xy)(xy)

    fig = plt.figure()
    fig, ax = using_mpl_scatter_density(
        fig, all_A_vs_B["logratio_Fitness_A"], all_A_vs_B["logratio_Fitness_B"]
    )
    ax.set_xlabel("$Fitness_{~Replicate~A~}$")
    ax.set_ylabel("$Fitness_{~Replicate~B~}$")
    corr, pval = stats.pearsonr(
        all_A_vs_B["logratio_Fitness_A"], all_A_vs_B["logratio_Fitness_B"]
    )
    ax.set_xlim([-5, 20])
    ax.set_ylim([-5, 20])
    ax.set_box_aspect(1)
    ax.annotate(
        f"Pearson's score: {round(corr,3)}\npval: {'{:.2e}'.format(pval)}",
        xy=[0, 0],
        xycoords="figure points",
    )
    plt.title(MTX)

    if not os.path.exists(f"{root}/figures/figure_1"):
        os.makedirs(f"{root}/figures/figure_1")

    fig.savefig(f"{root}/figures/figure_1/replicability_density_plot_{MTX}.eps", dpi=300)


def plot_number_strains_stat(logs, SINGLE, root):
    f = plt.figure(figsize=(SINGLE, SINGLE))
    nb_strain_stats = pd.concat(logs).reset_index(drop=True)
    nb_strain_stats["REP"] = [
        nb_strain_stats["Condition"][idx].split("_")[0] for idx in nb_strain_stats.index
    ]
    nb_strain_stats["MTX"] = [
        nb_strain_stats["Condition"][idx].split("_")[1] for idx in nb_strain_stats.index
    ]
    nb_strain_stats["DRUG"] = [
        nb_strain_stats["Condition"][idx].split(".")[0].split("_")[2]
        for idx in nb_strain_stats.index
    ]
    nb_strain_stats["DRUG"] = nb_strain_stats["DRUG"].astype(DRUG_ORDER)

    #### Stats

    print(len(nb_strain_stats))

    print(nb_strain_stats.describe())

    print(
        "Average number of strains per sample", np.mean(nb_strain_stats["nb_strains"])
    )

    print(
        len(nb_strain_stats[nb_strain_stats["nb_strains"] == 361])
        / len(nb_strain_stats)
    )

    print(
        len(
            nb_strain_stats[
                (nb_strain_stats["nb_strains"] >= 325)
                & (nb_strain_stats["nb_strains"] <= 360)
            ]
        )
        / len(nb_strain_stats)
    )

    sns.histplot(
        nb_strain_stats["nb_strains"],
        stat="percent",
        binwidth=10,
        binrange=(0, 361),
        color="grey",
    )
    plt.xticks(
        [0, 60, 120, 180, 240, 300, 360], ["0", "60", "120", "180", "240", "300", "360"]
    )
    plt.ylabel("Percentage (%)")
    plt.xlabel("Number of retrieved strains")
    f.savefig(f"{root}/figures/figure_1/strain_mapping_distribution.eps", dpi=300)


if __name__ == "__main__":
    curr_path = os.path.abspath(os.path.join(__file__, "../"))
    root_path = os.path.abspath(os.path.join(__file__, "../../.."))

    lineage_tracking_stats = pd.read_csv(f"{curr_path}/lineage_tracking_stats.csv")
    all_A_vs_B = pd.read_csv(f"{curr_path}/A_vs_B_annotated.csv")

    print(len(all_A_vs_B))

    #### FIGURE PANEL G
    # display_density_plot(all_A_vs_B, "", SINGLE, root_path)
    MTX = all_A_vs_B[all_A_vs_B["MTX_condition"].str.contains("_MTX")].reset_index(
        drop=True
    )
    display_density_plot(MTX, "MTX", SINGLE, root_path)
    # noMTX = all_A_vs_B[all_A_vs_B["MTX_condition"].str.contains("_noMTX")].reset_index(
    #     drop=True
    # )
    # display_density_plot(noMTX, "noMTX", SINGLE, root_path)

    #### FIGURE PANEL F
    log_folder = f"{root_path}/results/01_barcode_count/log"
    logs = [
        pd.read_csv(os.path.join(log_folder, log)) for log in os.listdir(log_folder)
    ]
    plot_number_strains_stat(logs, SINGLE, root_path)
