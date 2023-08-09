import numpy as np 
import pandas as pd 
import seaborn as sns
import matplotlib as mpl
import scipy.stats as stats
import matplotlib.pyplot as plt

from matplotlib_venn import venn2
from matplotlib_venn import venn3
from IPython.display import display

from venn import venn


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

CM = 1/2.54 # centimeters in inches
QUARTER = 4.5 * CM
THIRD = 6 * CM
HALF = 9 * CM
WHOLE = 18 * CM


#### Loading all SNPs annotations
all_ORF_SNPs = pd.read_csv('../data/genotype_information/snps_annotations_genome-version-3-64-1.txt')
print('Number of SNPs on ORF:', len(all_ORF_SNPs[all_ORF_SNPs['snps_class_up'] == 'ORF']))
ORF_with_SNPs = all_ORF_SNPs[all_ORF_SNPs['snps_class_up'] == 'ORF']['locus_id'].values

#### Loading human and yeast orthologs
human_homologs = pd.read_csv('../data/gwas_annotations/Yeast_vs_Human_Homologs.txt', sep='\t')
yeast_homologs_with_snps = set(human_homologs['Yeast_Gene_stable_ID']) & set(ORF_with_SNPs)
print(len(yeast_homologs_with_snps))

#### Loading piQTL results 
piQTLs = pd.read_csv('../results/05_piQTL_tables/significant_piQTLs_results_with_genome_annotations_without_MTX_peaks.csv')
ES_piQTLs = pd.read_csv('../results/05_piQTL_tables/ES_vs_NES/Metformin_ES_piQTLs.csv')
NES_piQTLs = pd.read_csv('../results/05_piQTL_tables/ES_vs_NES/Metformin_NES_piQTLs.csv')

#### Loading GWAS hits on diabetes type II
diabetes_typeII_gwas = pd.read_csv('../data/gwas_annotations/diabetes_type2-genes.csv')
human_diabetesII_QTLs = diabetes_typeII_gwas[~diabetes_typeII_gwas['Human_Gene_name'].str.contains('-')].drop_duplicates(('Human_Gene_name'))['Human_Gene_name']
yeast_diabetes_II_homologs = human_homologs[human_homologs['Human_Gene_name'].isin(human_diabetesII_QTLs)]['Yeast_Gene_stable_ID']
yeast_diabetes_II_homologs_with_snps = set(yeast_diabetes_II_homologs) & set(ORF_with_SNPs)

yeast_ES_piQTLs = [ gene for gene in np.unique(ES_piQTLs[(ES_piQTLs['snps_class_up'] == 'ORF')]['locus_id']) ]

yeast_NES_piQTLs = [ gene for gene in np.unique(NES_piQTLs[(NES_piQTLs['snps_class_up'] == 'ORF')]['locus_id']) ]
yeast_all_piQTLs = [ gene for gene in np.unique(piQTLs[(piQTLs['snps_class_up'] == 'ORF')]['locus_id']) ]

yeast_ES_piQTLs_with_homologs = set(yeast_ES_piQTLs) & set(yeast_homologs_with_snps)


f = plt.figure(figsize=(HALF, HALF))
out = venn3([set(yeast_homologs_with_snps), set(yeast_ES_piQTLs_with_homologs), set(yeast_diabetes_II_homologs_with_snps)], set_colors =('green','white', 'purple'), set_labels=['ORFs with annotated SNPs','Metformin \nORF piQTLs', 'Yeast ORF homologs of human genes \nassociated to Type II Diabetes'])
plt.show() 
plt.title('Drug: Metformin\n Indication:Type II Diabetes')
# f.savefig('../manuscript/figures/FIGURE_3/FIGURE_3A.pdf', dpi=300)
# f.savefig('../manuscript/figures/FIGURE_3/FIGURE_3A.eps', dpi=300)


N = 1290
n = 74
K = 89
k = 7
print(stats.hypergeom.pmf(k, N, K, n))

  
# creating data
data = [[7, 89], [74, 1290]]
  
# performing fishers exact test on the data
odd_ratio, p_value = stats.fisher_exact(data, alternative='two-sided')
print('odd ratio is : ' + str(odd_ratio))
print('p_value is : ' + str(p_value))

res = []
for locus_id in (set(yeast_ES_piQTLs) & set(yeast_diabetes_II_homologs)):
    print(locus_id)
    TMP = piQTLs[(piQTLs['locus_id'] == locus_id) & (piQTLs['Condition'].str.contains('Metformin'))]
    if len(TMP) != 0 :
        print(locus_id)
        res.append(TMP[['SNP', 'Condition', 'Chr', 'Pos', '-log_pval', 'locus_id', 'name', 'description']])


TABLE = pd.concat(res)
TABLE.to_csv('../results/06_GWAS_overlapping/human_GWAS_vs_piQTLs_metformin_hits.csv', index=False)
TABLE = TABLE.merge(human_homologs, left_on='locus_id', right_on='Yeast_Gene_stable_ID')
# TABLE.to_csv('../results/06_GWAS_overlapping/human_GWAS_vs_piQTLs_metformin_hits_with_human_homologs.csv', index=False)