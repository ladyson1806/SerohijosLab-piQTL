
import os 
import pandas as pd

from tqdm import tqdm 

#### Load MTX-specific piQTLs (piQTLs present in more than 155 conditions)
MTX_QTLs = pd.read_csv('../../results/05_piQTL_tables/MTX_specific_QTLs_with_genome_annotations.csv')
MTX_QTLs['SNP_ID'] = [ f'CHR{MTX_QTLs["Chr"][idx]}_{MTX_QTLs["SNP"][idx]}' for idx in MTX_QTLs.index ]
MTX_piQTLs = list(MTX_QTLs['SNP_ID'])

#### MTX-specific piQTLs identified from heatmap
MTX_piQTLs.append('CHR15_10435')
MTX_piQTLs.append('CHR15_10455')
MTX_piQTLs.append('CHR15_10468')


piQTLs_original = '../../data/piQTLs/piQTL_original/'
piQTLs_without_MTX = '../../data/piQTLs/piQTL_without_MTX/'

for piQTL_table in tqdm(os.listdir(piQTLs_original)):
    table = pd.read_csv(os.path.join(piQTLs_original, piQTL_table))
    table['SNP_ID'] = [ f'CHR{table["CHR"][idx]}_{table["SNP"][idx]}' for idx in table.index ]
    table_without_MTX_piQTLs = table[~table['SNP_ID'].isin(MTX_piQTLs)] 
    table_without_MTX_piQTLs.to_csv(f'{piQTLs_without_MTX}{piQTL_table}', index=False)




