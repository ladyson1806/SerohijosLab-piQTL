import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import cv2, os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import mpl_scatter_density
import scipy.stats as stats
import matplotlib.pyplot as plt


from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
from pandas.api.types import CategoricalDtype
from matplotlib.colors import LinearSegmentedColormap


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


if __name__ == "__main__":
    all_A_vs_B = pd.read_csv("./pipeline/A_vs_B_annotated_replicates.csv")

    print(len(all_A_vs_B))

    #### FIGURE EXT FIG 3A
    # print('Best replicate')
    # print(lineage_tracking_stats[lineage_tracking_stats['Correlation'] == lineage_tracking_stats['Correlation'].max()])
    # # fitness_comparison('HUG1_PGM2', '5-FC', DOUBLE)

    #### FIGURE EXT FIG 3B
    # display_density_plot(all_A_vs_B, '', SINGLE)
    MTX = all_A_vs_B[all_A_vs_B["MTX_condition"].str.contains("_MTX")].reset_index(
        drop=True
    )
    # display_density_plot(MTX, 'MTX', SINGLE)
    noMTX = all_A_vs_B[all_A_vs_B["MTX_condition"].str.contains("_noMTX")].reset_index(
        drop=True
    )
    # display_density_plot(noMTX, 'noMTX', SINGLE)

    i = 0
    title = ["", "MTX", "noMTX"]
    for TABLE in [all_A_vs_B, MTX, noMTX]:
        REP1 = TABLE[TABLE["strain_id"] == "17"]
        REP2 = TABLE[TABLE["strain_id"] == "17_2"]
        REP_17 = REP1.merge(REP2, on="MTX_condition", suffixes=["_REP1", "_REP2"])
        REP_17["REP_ID"] = "17"

        REP1 = TABLE[TABLE["strain_id"] == "40"]
        REP2 = TABLE[TABLE["strain_id"] == "40_2"]
        REP_40 = REP1.merge(REP2, on="MTX_condition", suffixes=["_REP1", "_REP2"])
        REP_40["REP_ID"] = "40"

        REP1 = TABLE[TABLE["strain_id"] == "180"]
        REP2 = TABLE[TABLE["strain_id"] == "180_2"]
        REP_180 = REP1.merge(REP2, on="MTX_condition", suffixes=["_REP1", "_REP2"])
        REP_180["REP_ID"] = "180"

        REP1 = TABLE[TABLE["strain_id"] == "43_ref1"]
        REP2 = TABLE[TABLE["strain_id"] == "43_ref2"]
        REP_43_REF = REP1.merge(REP2, on="MTX_condition", suffixes=["_REP1", "_REP2"])
        REP_43_REF["REP_ID"] = "43_REF"

        REP1 = TABLE[TABLE["strain_id"] == "599_ref1"]
        REP2 = TABLE[TABLE["strain_id"] == "599_ref2"]
        REP_599_REF = REP1.merge(REP2, on="MTX_condition", suffixes=["_REP1", "_REP2"])
        REP_599_REF["REP_ID"] = "599_REF"

        REPS = pd.concat([REP_17, REP_40, REP_180, REP_43_REF, REP_599_REF])

        f = plt.figure(figsize=(12 * CM, 12 * CM))
        sns.scatterplot(
            x=REPS["logratio_Fitness_A_REP1"],
            y=REPS["logratio_Fitness_A_REP2"],
            hue=REPS["REP_ID"],
            marker=".",
            linewidth=0,
            alpha=0.7,
            legend=False,
            palette=["blue", "orange", "green", "red", "purple"],
        )
        sns.scatterplot(
            x=REPS["logratio_Fitness_B_REP1"],
            y=REPS["logratio_Fitness_B_REP2"],
            hue=REPS["REP_ID"],
            marker=".",
            linewidth=0,
            alpha=0.7,
            legend=False,
            palette=["blue", "orange", "green", "red", "purple"],
        )

        REPS_1 = []
        for fitness in REPS["logratio_Fitness_A_REP1"]:
            REPS_1.append(fitness)
        for fitness in REPS["logratio_Fitness_B_REP1"]:
            REPS_1.append(fitness)

        REPS_2 = []
        for fitness in REPS["logratio_Fitness_A_REP2"]:
            REPS_2.append(fitness)
        for fitness in REPS["logratio_Fitness_B_REP2"]:
            REPS_2.append(fitness)

        corr, pval = stats.pearsonr(REPS_1, REPS_2)
        plt.xlabel("$Fitness_{~Duplicate~1}$")
        plt.ylabel("$Fitness_{~Duplicate~2}$")
        plt.xlim([-15, 8])
        plt.ylim([-15, 8])
        plt.annotate(
            f"Pearson's corr: {round(corr,3)}\npval: {'{:.2e}'.format(pval)}",
            xy=[0, 270],
            xycoords="axes points",
        )

        plt.title(title[i])

        custom = [
            Line2D([], [], marker=".", color="blue", linestyle="None"),
            Line2D([], [], marker=".", color="orange", linestyle="None"),
            Line2D([], [], marker=".", color="green", linestyle="None"),
            Line2D([], [], marker=".", color="red", linestyle="None"),
            Line2D([], [], marker=".", color="purple", linestyle="None"),
        ]

        plt.legend(
            custom,
            ["17", "40", "180", "43_REF", "599_REF"],
            fontsize=7,
            loc="upper left",
        )

        f.savefig(
            f"../manuscript/figures/EXT_FIGURE_3/duplicate_density_plot_{title[i]}.png",  # ToDo: modify output path
            dpi=300,
        )
        f.savefig(
            f"../manuscript/figures/EXT_FIGURE_3/duplicate_density_plot_{title[i]}.eps",  # ToDo: modify output path
            dpi=300,
        )
        i += 1
