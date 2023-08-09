import matplotlib.pyplot as plt 
import scipy.stats as stats
import seaborn as sns
import plotnine as p9
import pandas as pd
import numpy as np
import allel, re

from pandas.api.types import CategoricalDtype
from scipy.spatial.distance import squareform
from scipy.io import loadmat
from tqdm import tqdm

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks")
sns.set(font="Helvetica")

plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams.update({'font.size': 7})


chr_order = CategoricalDtype(
    ['CHR_1', 
    'CHR_2', 
    'CHR_3', 
    'CHR_4', 
    'CHR_5', 
    'CHR_6', 
    'CHR_7', 
    'CHR_8', 
    'CHR_9', 
    'CHR_10', 
    'CHR_11', 
    'CHR_12', 
    'CHR_13', 
    'CHR_14', 
    'CHR_15', 
    'CHR_16', 
    'CHR_MT'], 
    ordered=True
)

def get_chr_limits(snp_position):
    #### Delimitations of chromosomes based on SNPs index
    chr_limits = []
    for chr in np.unique(snp_position['chrom'].values):
        TMP = snp_position[snp_position['chrom'] == chr]
        chr_limits.append([chr, int(TMP['SNP'].describe()['min']), int(TMP['SNP'].describe()['max']), int(TMP['position'].describe()['min']), int(TMP['position'].describe()['max'])])
    chr_limit_table = pd.DataFrame(chr_limits, columns=['CHR', 'SNP_IDX_FIRST', 'SNP_IDX_LAST', 'SNP_POS_FIRST', 'SNP_POS_LAST']).sort_values('SNP_IDX_FIRST')
    chr_limit_table['CHR'] = chr_limit_table['CHR'].astype(chr_order)
    print(chr_limit_table)
    return chr_limit_table

def create_genotype_map(genotype_matrix):
#### Remove snp_id column
    cols = genotype_matrix.columns
    genotype_matrix = genotype_matrix[cols[1:]]
    print(genotype_matrix.shape)
    f, ax = plt.subplots(nrows=1, ncols=1)
    ax.grid(False)
    plt.imshow(genotype_matrix.T.values, vmin=-1, vmax=1, extent=[1, 12054, 357, 1], cmap='viridis', interpolation='None', aspect='auto')
    ### Change SNP index to chromosome label
    # ax.set_xticks([1, 182, 799, 1619, 2610, 3456, 4210, 5248, 6219, 6924, 7290, 8805, 9629, 9894, 10879, 11472, 12020])
    # chr_order = list(chr_limits.sort_values('CHR')['CHR'])
    # ax.set_xticklabels(chr_order, rotation=90, fontsize=7)
    plt.colorbar()

    f.savefig('../manuscript/figures/EXT_FIGURE_1/FIGURE_1A.png', dpi=300)
    f.savefig('../manuscript/figures/EXT_FIGURE_1/FIGURE_1A.eps', dpi=300)


def format_genotype(x):
    if x == -1 :
        return 0
    elif x == 1 :
        return 2
    elif x == 0 :
        return -1

def linkage_disequilibrium(genotype_matrix):
#### LD with allel library
    # cols = genotype_matrix.columns
    # gn = genotype_matrix[cols[1:]].T.apply(np.vectorize(format_genotype)).T.values
    # r = allel.rogers_huff_r(gn)
    # r_table = pd.DataFrame(squareform(r ** 2))
    # r_table.columns = [ col+1 for col in r_table.columns ]
    # r_table.index = range(1, len(r_table) + 1)
    genotype_matrix = genotype_matrix.replace(0, "NA")
    print(genotype_matrix)

    f = open('../data/genotype_information/piQTL_genotype_matrix_for_LD_dec2022.txt', "w")
    for idx in tqdm(genotype_matrix.index) : 
        for col in genotype_matrix.columns[1:] :
            line = '\t'.join([str(genotype_matrix['snp_id'][idx]), col, str(genotype_matrix.loc[idx, col]), str(1)])
            f.write(line+'\n')
    f.close()
    
    # f, ax = plt.subplots(nrows=1, ncols=1)
    # plt.imshow(r_table.values, vmin=0, vmax=1, extent=[1, 12054, 12054, 1], cmap='binary', interpolation='none', aspect='equal')
    # ### Change SNP index to chromosome label
    # ax.set_xticks([1, 182, 799, 1619, 2610, 3456, 4210, 5248, 6219, 6924, 7290, 8805, 9629, 9894, 10879, 11472, 12020])
    # ax.set_yticks([1, 182, 799, 1619, 2610, 3456, 4210, 5248, 6219, 6924, 7290, 8805, 9629, 9894, 10879, 11472, 12020])
    # chr_order = list(chr_limits.sort_values('CHR')['CHR'])
    # ax.set_xticklabels(chr_order, rotation=90, fontsize=10)
    # ax.set_yticklabels(chr_order, fontsize=10)
    # plt.colorbar()
    # plt.show()
    # f.savefig('../manuscript/figures/FIGURE_1/FIGURE_1B.png', dpi=300)
    # f.savefig('../manuscript/figures/FIGURE_1/FIGURE_1B.svg', dpi=300)

    # return r_table

# def plot_LD_per_chromosome(r_table, chr):
#     pass 

# def get_all_snp_distance(snp_position):
#     all_snp_distance = []
#     for chr in np.unique(snp_position['chrom'].values):
#         TMP = snp_position[snp_position['chrom'] == chr].sort_values('position')
#         chr_positions = list(TMP['position'].values)
#         for i in range(len(chr_positions)-1) :
#             snp_distance = chr_positions[i+1] - chr_positions[i]
#             all_snp_distance.append([chr, snp_distance])
#     snp_distance_table = pd.DataFrame(all_snp_distance, columns=['CHR', 'distance'])
#     return snp_distance_table


if __name__ == "__main__":
    genotype_matrix = pd.read_csv('../data/genotype_information/piQTL_genotype_matrix_dec2022.txt')
    snp_position = pd.read_csv('../data/genotype_information/snps_annotations_genome-version-3-64-1.txt').rename(columns={'snp_id':'SNP'})
    print(genotype_matrix.shape)
    chr_limits = get_chr_limits(snp_position)
    print("Save panel A for Figure 1")
    create_genotype_map(genotype_matrix)
    print("Input for generating panel B for Figure 1 with Rscript calculate_LD.R ")
    # linkage_disequilibrium(genotype_matrix)

    

    # #### Distribution of SNPs distance (in bp) per chromosome
    # snp_distance_table = get_all_snp_distance(snp_position)
    # print(snp_distance_table.describe(percentiles=[0.1,0.5,0.9]))
    # CHR_order = [ f'CHR{i+1}' for i in range(16)]
    # CHR_order.append('CHR_MT')
    # no_extreme_outliers = snp_distance_table[snp_distance_table['distance'] < 3000]
    # boxplot_distance = (
    #     p9.ggplot(no_extreme_outliers, p9.aes(x='CHR', y='distance'))
    #     + p9.geom_boxplot()
    #     + p9.scale_x_discrete(limits=CHR_order)
    #     + p9.labs(y='Distance (in bp) separating two SNPs', x='Chromosome')
    #     + p9.theme_classic()
    #     + p9.theme(axis_text_x = p9.element_text(angle=90))
    # )
    # boxplot_distance.save('./FIGURES/SNPs_distance_distribution_per_chrom.svg', dpi=200)
    # boxplot_distance.save('./FIGURES/SNPs_distance_distribution_per_chrom.png', dpi=200)
    
    # #### Distribution of SNPs distance (in bp) without chromosome distinction
    # g = plt.figure(figsize=(6,6))
    # sns.displot(x=no_extreme_outliers['distance'])
    # plt.ylabel('Count')
    # plt.xlabel('Distance (in bp) separating two SNPs')
    # plt.show()
