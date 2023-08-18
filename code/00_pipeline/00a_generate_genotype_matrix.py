import re
from os import path

import numpy as np
import pandas as pd
import roman
from pandas.api.types import CategoricalDtype
# from scipy.io import loadmat

chr_order = CategoricalDtype(
    [
        'CHR_1',
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
        'CHR_MT'
    ],
    ordered=True
)


def get_position_value(x):
    try:
        pos = re.findall(r'[0]+([0-9]+)', x)[0]
        return int(pos)
    except:
        return int(x)


def format_snp_id(x):
    x = x.replace('chr', 'snp')
    return x


def format_roman_chrom(x):
    if 'Mito' not in x :
        rom_nb = x
        nb = roman.fromRoman(rom_nb)
        x = 'CHR_'+ str(nb)
        return x
    return 'CHR_MT'


def format_gene_loc(gene_loc, root):
    # gene_loc = gene_loc[['Gene', 'Chromosome_Name', 'start', 'end']]
    # gene_loc = gene_loc.rename(columns={'Gene':'gene_id', 'Chromosome_Name':'chrom'})
    gene_loc = gene_loc[gene_loc['genome_annotations'].str.contains('ORF')]
    gene_loc[['locus_id', 'name', 'sgd_id', 'chrom', 'start', 'end']].to_csv(f'{root}/data/genotype_information/yeast_ORFs_dec2022.txt', index=False)


def format_snps_loc(variant_pos, root):
    # chrom_table = pd.read_csv('../../data/genotype_information/genes_by_strain/S288C_chromosomes.csv', sep=',').rename(columns={'chromosomeID':'ncbi_id_ref', 'chromosome_label':'chrom'})
    # chrom_table['ncbi_id_ref'] = [ chrom_table['ncbi_id_ref'][i].split('.')[0] for i in chrom_table.index ]
    # chrom_table = chrom_table[['ncbi_id_ref', 'chrom']]

    # var_pos_list = []
    # for i in range(len(variant_pos['variantPos'])):
    #     variant = variant_pos['variantPos'][i][0][0]
    #     tmp = []
    #     for var in variant.split(':'):
    #         if '|' in var :
    #             var = var.split('|')[1]
    #             tmp.append(var)
    #         else:
    #             tmp.append(var)
    #     var_pos_list.append(tmp)

    # #### Create snp table
    # var_pos_table = pd.DataFrame(var_pos_list, columns=['ncbi_id_ref', 'position', 'REF', 'ALT', 'GT_REF', 'GT_ALT'])
    # var_pos_table = var_pos_table.merge(chrom_table, on='ncbi_id_ref')
    # #### Formatting chromosome, snp_id, positions columns
    # var_pos_table['genetic_distance'] = 0
    # var_pos_table['position'] = var_pos_table['position'].apply(get_position_value)
    # var_pos_table['chrom'] = var_pos_table['chrom'].astype(chr_order)
    # var_pos_table = var_pos_table.sort_values(['chrom', 'position'])
    # var_pos_table['snp_id'] = [ x+1 for x in range(len(var_pos_table)) ]
    # var_pos_table[['snp_id', 'chrom', 'position', 'REF', 'ALT']].to_csv('../../data/genotype_information/yeast_snps_loc_test.txt', index=False)
    # return var_pos_table

    variant_pos['snp_id'] = variant_pos.index + 1
    variant_pos['chrom'] = variant_pos['chrom'].apply(format_roman_chrom)
    variant_pos[['snp_id', 'chrom', 'position', 'REF', 'ALT']].to_csv(f'{root}/data/genotype_information/yeast_snps_loc_dec2022.txt', index=False)



def format_genotypic_matrix(genotype_matrix, barcoded_strains, snps_loc, root):
    ###### Previous version with loadmat() from scipy
    # genotype_table = pd.DataFrame(genotype_matrix['phasedGenotype'])
    # genotype_table = genotype_table.astype(int)

    # #### Adding strain_ids
    # genotype_table['strain_id'] = [ x+1 for x in range(len(genotype_matrix['phasedGenotype']))]

    # #### Transposing table and removing the line used for adding strain_ids as columns names
    # genotype_table = genotype_table.T.reset_index(drop=True)
    # genotype_table.columns = genotype_table.iloc[12054]
    # cols = list(np.unique(barcoded_strains))
    # genotype_table = genotype_table[cols]
    # genotype_table = genotype_table.drop(12054)
    # genotype_table.columns = genotype_table.columns.astype(str)

    # #### Adding duplicated strains
    # genotype_table['17_2'] = genotype_table['17']
    # genotype_table['40_2'] = genotype_table['40']
    # genotype_table['180_2'] = genotype_table['180']
    # # genotype_table['43_ref1'] = genotype_table['43']
    # # genotype_table['599_ref1'] = genotype_table['599']
    # # genotype_table['43_ref2'] = genotype_table['43']
    # # genotype_table['599_ref2'] = genotype_table['599']

    # #### Addind snp_ids column and put at the first position
    # genotype_table['snp_id'] = snps_loc['snp_id']

    #############
    #### Adding strain_ids
    genotype_matrix['strain_id'] = [ x+1 for x in genotype_matrix.index ]
    #### Transposing table and removing the line used for adding strain_ids as columns names
    genotype_table = genotype_matrix.T.reset_index(drop=True)
    genotype_table.columns = genotype_table.iloc[12054]
    cols = list(np.unique(barcoded_strains))
    genotype_table = genotype_table[cols]
    genotype_table = genotype_table.drop(12054)
    genotype_table.columns = genotype_table.columns.astype(str)

    #### Adding duplicated strains
    genotype_table['17_2'] = genotype_table['17']
    genotype_table['40_2'] = genotype_table['40']
    genotype_table['180_2'] = genotype_table['180']

    genotype_table =  genotype_table.reset_index().rename(columns={'index':'snp_id'})
    genotype_table['snp_id'] = genotype_table['snp_id'] + 1
    genotype_table.to_csv(f'{root}/data/genotype_information/piQTL_genotype_matrix_dec2022.txt', index=False)


if __name__ == "__main__":
    root_path = path.abspath(path.join(__file__, "../../.."))

    ### Required inputs
    genotype_matrix = pd.read_csv(f'{root_path}/data/genotype_information/she_jarosz_2018/phasedGenotype_from_mat.csv', header=None)
    variant_position_matrix = pd.read_csv(f'{root_path}/data/genotype_information/she_jarosz_2018/variantPos_from_vcf.csv', sep='\t', names=['chrom', 'position', 'REF', 'ALT'])
    gene_loc = pd.read_csv(f'{root_path}/data/genome_annotations/sgd_database/orf_coding_R64-3-1.csv')

    #### List of barcoded strains
    july_barcoded_strain_table = pd.read_csv(f'{root_path}/data/pipeline/barcode_collection.csv')
    barcoded_strains = list(july_barcoded_strain_table['strain_number'])

    #### Save gene (CDS) information table
    gene_loc = format_gene_loc(gene_loc, root_path)

    #### Save snps information table
    snps_loc = format_snps_loc(variant_position_matrix, root_path)

    #### Genotype of piQTL subset from She and Jarosz study
    genotype_table = format_genotypic_matrix(genotype_matrix, barcoded_strains, snps_loc, root_path)
