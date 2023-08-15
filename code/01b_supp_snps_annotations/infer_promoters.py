import numpy as np
import pandas as pd
import plotnine as p9


def build_promoter_table(minus_one, plus_one):
    nucleosome_positions = minus_one.merge(
        plus_one, on=["chrom", "gene", "orientation"], suffixes=("_-1", "_+1")
    )[
        [
            "chrom",
            "gene",
            "stable/fragile",
            "start_-1",
            "end_-1",
            "start_+1",
            "end_+1",
            "orientation",
        ]
    ]

    gene_plus_oriented = nucleosome_positions[
        nucleosome_positions["orientation"] == "+"
    ].reset_index(drop=True)
    gene_plus_oriented["NFR_start"] = gene_plus_oriented["end_-1"] + 1
    gene_plus_oriented["NFR_end"] = gene_plus_oriented["start_+1"] - 1
    gene_plus_oriented["NFR_length"] = (
        gene_plus_oriented["NFR_end"] - gene_plus_oriented["NFR_start"]
    )
    gene_plus_oriented["promoter_id"] = (
        gene_plus_oriented["gene"] + "_" + gene_plus_oriented["stable/fragile"]
    )

    gene_minus_oriented = nucleosome_positions[
        nucleosome_positions["orientation"] == "-"
    ].reset_index(drop=True)
    gene_minus_oriented["NFR_start"] = gene_minus_oriented["start_+1"] - 1
    gene_minus_oriented["NFR_end"] = gene_minus_oriented["end_-1"] + 1
    gene_minus_oriented["NFR_length"] = (
        gene_minus_oriented["NFR_end"] - gene_minus_oriented["NFR_start"]
    )
    gene_minus_oriented["promoter_id"] = (
        gene_minus_oriented["gene"] + "_" + gene_minus_oriented["stable/fragile"]
    )

    format_chromosome = dict()
    for i in range(16):
        format_chromosome[np.unique(nucleosome_positions["chrom"])[i]] = i + 1

    NFR_table = pd.concat([gene_plus_oriented, gene_minus_oriented]).reset_index(
        drop=True
    )
    NFR_table["chrom"] = [
        format_chromosome[NFR_table["chrom"][i]] for i in NFR_table.index
    ]
    print(
        f'Removing promoters with NFR length > 1000: {len(NFR_table[NFR_table["NFR_length"] > 1000])}'
    )
    NFR_table = NFR_table[NFR_table["NFR_length"] <= 1000].reset_index(drop=True)
    return NFR_table


if __name__ == "__main__":
    minus_one = pd.read_csv(
        "./data/nucleosome_position/minus_one_nucleosome_positions.csv"
    ).rename(columns={"chromosome": "chrom", "gene orientation": "orientation"})
    plus_one = pd.read_csv(
        "./data/nucleosome_position/plus_one_nucleosome_positions.csv"
    ).rename(columns={"chromosome": "chrom"})
    promoter_position = build_promoter_table(minus_one, plus_one)
    # promoter_position.to_csv('./data/yeast_promoter_loc.txt', index=False)
    print("Promoter location table generated")

    print("About the yeast promoters")

    print(np.unique(promoter_position["orientation"], return_counts=True))

    NFR_COUNTPLOT = (
        p9.ggplot(data=promoter_position)
        + p9.geom_bar(p9.aes(x="stable/fragile"))
        + p9.geom_text(
            p9.aes(x="stable/fragile", label=p9.after_stat("count")),
            stat="count",
            nudge_y=0.125,
            va="bottom",
        )
        + p9.labs(x="Promoter type")
        + p9.theme_classic()
    )
    NFR_COUNTPLOT.save("./figures/NFR_counts.png", dpi=150)

    NFR_BOXPLOT = (
        p9.ggplot(data=promoter_position)
        + p9.geom_violin(p9.aes(y="NFR_length", x="stable/fragile"))
        + p9.labs(x="Promoter type", y="Promoter length")
        + p9.theme_classic()
    )
    NFR_BOXPLOT.save("./figures/NFR_length_violinplot.png", dpi=150)
