import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

import os
from functools import reduce

import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm


def reads_per_ppi(barcode_count_folder, PPI_list, results_folder):
    tables = [ pd.read_csv(os.path.join(barcode_count_folder, table)) for table in os.listdir(barcode_count_folder) if '.csv' in table ]
    ALL = reduce(lambda left,right: pd.merge(left,right,on=['strain_id'], how='outer'), tables).replace(np.nan, 0)

    for PPI in sorted(PPI_list):
        cols = sorted([ col for col in ALL.columns if PPI in col ])
        cols.insert(0, 'strain_id')
        ALL[cols].to_csv(f'{results_folder}/01_barcode_count/per_PPI/{PPI}_read_number.csv', index=False)


def logratio_fitness(PPI, DRUG, data_folder):
    ALL_MTX = []
    for MTX in ['noMTX', 'MTX']:
        ALL_REP = []
        for REP in ['A', 'B']:
            PPI_table = pd.read_csv(f'{data_folder}/results/01_barcode_count/per_PPI/{PPI}_read_number.csv')
            R = PPI_table[['strain_id', f'T0_noMTX_noDrug.{PPI}', f'{REP}_{MTX}_{DRUG}.{PPI}']]
            R[f'{REP}_{MTX}_{DRUG}.{PPI}'] = R[f'{REP}_{MTX}_{DRUG}.{PPI}'] + 1
            R[f'T0.{PPI}'] = R[f'T0_noMTX_noDrug.{PPI}'] + 1

            R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] = (R[f'{REP}_{MTX}_{DRUG}.{PPI}'] / sum(R[f'{REP}_{MTX}_{DRUG}.{PPI}'])) * 10**6
            R[f'T0.{PPI}_RPM'] = (R[f'T0.{PPI}'] / sum(R[f'T0.{PPI}'])) * 10**6
            R[f'{PPI}_{MTX}_{DRUG}_logratio_Fitness'] = np.log2(R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] / R[f'T0.{PPI}_RPM'])
            RPM_cols = [ col for col in R.columns if 'RPM' in col ]
            RPM_cols.insert(0, 'strain_id')
            RPM_cols.insert(-1, f'{PPI}_{MTX}_{DRUG}_logratio_Fitness')
            #  display(R[RPM_cols])
            ALL_REP.append(R[RPM_cols])

        #### Average Fitness from mean of logratio_A and logratio_B
        REP_MERGED = pd.merge(ALL_REP[0], ALL_REP[1], on=['strain_id', f'T0.{PPI}_RPM'], suffixes=['_A', '_B'])
        REP_MERGED[f'{PPI}_{MTX}_{DRUG}_avg_logratio_Fitness'] = (REP_MERGED[f'{PPI}_{MTX}_{DRUG}_logratio_Fitness_A'] + REP_MERGED[f'{PPI}_{MTX}_{DRUG}_logratio_Fitness_B']) / 2
        ALL_MTX.append(REP_MERGED[['strain_id', f'{PPI}_{MTX}_{DRUG}_avg_logratio_Fitness']])
    MTX_MERGED = pd.merge(ALL_MTX[0],  ALL_MTX[1], on='strain_id')
    return MTX_MERGED


def run_logratio_fitness_per_ppi(PPI_list, results_folder, data_folder):
    for PPI in tqdm(PPI_list) :
        ALL_DRUGS = []
        for DRUG in ['noDrug', '5-FC', 'Fluconazole', 'Metformin', 'Trifluoperazine'] :
            REP_MERGED = logratio_fitness(PPI, DRUG, data_folder)
            ALL_DRUGS.append(REP_MERGED)
        FINAL_TABLE = reduce(lambda  left,right: pd.merge(left,right,on=['strain_id'], how='outer'), ALL_DRUGS)

        logratio_cols = [ col for col in FINAL_TABLE.columns if '_logratio_Fitness' in col ]
        logratio_cols.insert(0, 'strain_id')
        # display(FINAL_TABLE[logratio_cols])
        output_file = os.path.join(results_folder, f'03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv')
        FINAL_TABLE[logratio_cols].to_csv(output_file, index=False)


def ppi_quantification_before_normalization(PPI_list, results_folder):
    PPI_table = [ pd.read_csv(f'{results_folder}/03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv') for PPI in PPI_list ]
    all_logratio_fitness = reduce(lambda  left,right: pd.merge(left,right,on=['strain_id'], how='outer'), PPI_table)
    all_logratio_fitness = all_logratio_fitness[~all_logratio_fitness['strain_id'].str.contains('_ref')].rename(columns={'strain_id':'TAXA'})

    logratio_output_file = os.path.join(results_folder, '03_ppi_estimation/logratio/all_PPI_logratio_fitness.csv')
    all_logratio_fitness.to_csv(logratio_output_file, index=False)

    PPI_metadata = []
    for cdt in all_logratio_fitness.columns:
        if cdt !='TAXA':
            PPI = f"{cdt.split('_')[0]}_{cdt.split('_')[1]}"
            MTX = cdt.split('_')[2]
            DRUG = cdt.split('_')[3]

            PPI_metadata.append([cdt, PPI, DRUG, MTX])

    PPI_metadata_table = pd.DataFrame(PPI_metadata, columns=['label', 'PPI', 'Drug', 'MTX'])
    metadata_output_file = os.path.join(results_folder, '03_ppi_estimation/logratio/all_PPI_logratio_metadata.csv')
    PPI_metadata_table.to_csv(metadata_output_file, index=False)


def ppi_quantification_after_normalization(PPI_list, results_folder, genotype):
    ### Substraction of mean values from ref strains
    custom_cols = ['strain_id', 'noMTX_noDrug_avg_logratio_Fitness', 'MTX_noDrug_avg_logratio_Fitness', 'noMTX_5-FC_avg_logratio_Fitness', 'MTX_5-FC_avg_logratio_Fitness', 'noMTX_Fluconazole_avg_logratio_Fitness', 'MTX_Fluconazole_avg_logratio_Fitness',
'noMTX_Metformin_avg_logratio_Fitness', 'MTX_Metformin_avg_logratio_Fitness', 'noMTX_Trifluoperazine_avg_logratio_Fitness', 'MTX_Trifluoperazine_avg_logratio_Fitness']

    PPI_table = [ pd.read_csv(f'{results_folder}/03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv', header=None, names=custom_cols).drop(0).reset_index(drop=True) for PPI in PPI_list ]

    noPPI_reference_table = PPI_table[-1]
    cols = list(noPPI_reference_table.columns)
    cols.remove('strain_id')
    for col in cols:
        noPPI_reference_table[col] = noPPI_reference_table[col].astype(float)
    ALL_PPI = []
    for i in range(len(PPI_table)):
        PPI = PPI_list[i]
        cols = list(PPI_table[i].columns)
        cols.remove('strain_id')
        for col in cols:
            PPI_table[i][col] = PPI_table[i][col].astype(float)

        ref_mean = pd.DataFrame(PPI_table[i][PPI_table[i]['strain_id'].str.contains('_ref')].mean()).T
        for col in PPI_table[i].columns :
            if 'strain_id' not in col :
                PPI_table[i][f'{col}_minus_ref'] = PPI_table[i][col] - ref_mean[col].values[0]

        custom_cols_ref = [ col for col in PPI_table[i].columns if '_ref' in col ]
        custom_cols_ref.insert(0, 'strain_id')

        delta_reference = PPI_table[i][custom_cols_ref]
        final_custom_cols = [f'{PPI}_{col}' for col in delta_reference.columns[1:] ]
        final_custom_cols.insert(0, 'strain_id')
        delta_reference.columns = final_custom_cols
        ALL_PPI.append(delta_reference)

    PPI_FINAL_TABLE_normalized = reduce(lambda  left,right: pd.merge(left,right,on=['strain_id'], how='outer'), ALL_PPI)
    PPI_FINAL_TABLE_normalized = PPI_FINAL_TABLE_normalized[~PPI_FINAL_TABLE_normalized['strain_id'].str.contains('_ref')].rename(columns={'strain_id':'TAXA'})

    ### eQTL_matrix
    PPI_FINAL_TABLE_normalized.rename(columns={'TAXA':'Condition'}).set_index('Condition').T[genotype.columns[1:]].to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_for_eQTL_matrix.csv')

    #### rMVP
    PPI_FINAL_TABLE_normalized.rename(columns={'TAXA':'Condition'}).set_index('Condition').T[genotype.columns[1:]].T.reset_index().rename(columns={"Condition":"TAXA"}).to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_for_rMVP.csv', index=False)

    PPI_metadata = []

    for cdt in PPI_FINAL_TABLE_normalized.columns:
        if cdt !='TAXA':
            PPI = f"{cdt.split('_')[0]}_{cdt.split('_')[1]}"
            MTX = cdt.split('_')[2]
            DRUG = cdt.split('_')[3]
            PPI_metadata.append([cdt, PPI, DRUG, MTX])

    PPI_metadata_table = pd.DataFrame(PPI_metadata, columns=['label', 'PPI', 'Drug', 'MTX'])
    PPI_metadata_table.to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_metadata.csv', index=False)

def ppi_quantification_after_normalization_v2(PPI_list, results_folder, genotype):
    custom_cols = ['strain_id', 'noMTX_noDrug_avg_logratio_Fitness', 'MTX_noDrug_avg_logratio_Fitness', 'noMTX_5-FC_avg_logratio_Fitness', 'MTX_5-FC_avg_logratio_Fitness', 'noMTX_Fluconazole_avg_logratio_Fitness', 'MTX_Fluconazole_avg_logratio_Fitness',
    'noMTX_Metformin_avg_logratio_Fitness', 'MTX_Metformin_avg_logratio_Fitness', 'noMTX_Trifluoperazine_avg_logratio_Fitness', 'MTX_Trifluoperazine_avg_logratio_Fitness']

    PPI_table = [ pd.read_csv(f'{results_folder}/03_ppi_estimation/logratio/before_downsampling/{PPI}_logratio_fitness.csv', header=None, names=custom_cols).drop(0).reset_index(drop=True) for PPI in PPI_list ]

    noPPI_reference_table = PPI_table[-1]
    cols = list(noPPI_reference_table.columns)
    cols.remove('strain_id')
    for col in cols:
        noPPI_reference_table[col] = noPPI_reference_table[col].astype(float)

    ALL_PPI = []
    for i in range(len(PPI_table)-1):
        PPI = PPI_list[i]
        # print(PPI)
        cols = list(PPI_table[i].columns)
        cols.remove('strain_id')
        for col in cols:
            PPI_table[i][col] = PPI_table[i][col].astype(float)
        # display(PPI_table[i])
        ref_mean = pd.DataFrame(PPI_table[i][PPI_table[i]['strain_id'].str.contains('_ref')].mean()).T
        # display(ref_mean)
        for col in PPI_table[i].columns :
            if 'strain_id' not in col :
                PPI_table[i][f'{col}'] = PPI_table[i][col] - ref_mean[col].values[0]
        # display(delta_reference)
        delta_reference = PPI_table[i].set_index('strain_id').subtract(noPPI_reference_table.set_index('strain_id')).reset_index()
        final_custom_cols = [f'{PPI}_{col}' for col in delta_reference.columns[1:] ]
        final_custom_cols.insert(0, 'strain_id')
        delta_reference.columns = final_custom_cols
        ALL_PPI.append(delta_reference)

    PPI_FINAL_TABLE_normalized = reduce(lambda  left,right: pd.merge(left,right,on=['strain_id'], how='outer'), ALL_PPI)
    PPI_FINAL_TABLE_normalized = PPI_FINAL_TABLE_normalized[~PPI_FINAL_TABLE_normalized['strain_id'].str.contains('_ref')].rename(columns={'strain_id':'TAXA'})
    # display(PPI_FINAL_TABLE_normalized)
    PPI_FINAL_TABLE_normalized.rename(columns={'TAXA':'Condition'}).set_index('Condition').T[genotype.columns[1:]].to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_minus_noPPI-labelled_for_eQTL_matrix.csv')

    #### rMVP
    PPI_FINAL_TABLE_normalized.rename(columns={'TAXA':'Condition'}).set_index('Condition').T[genotype.columns[1:]].T.reset_index().rename(columns={"Condition":"TAXA"}).to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_noPPI-labelled_for_rMVP.csv', index=False)

    PPI_metadata = []

    for cdt in PPI_FINAL_TABLE_normalized.columns:
        if cdt !='TAXA':
            PPI = f"{cdt.split('_')[0]}_{cdt.split('_')[1]}"
            MTX = cdt.split('_')[2]
            DRUG = cdt.split('_')[3]

            PPI_metadata.append([cdt, PPI, DRUG, MTX])

    PPI_metadata_table = pd.DataFrame(PPI_metadata, columns=['label', 'PPI', 'Drug', 'MTX'])
    PPI_metadata_table.to_csv(f'{results_folder}/03_ppi_estimation/logratio/all_PPI_logratio_minus_ref_noPPI-labelled_metadata.csv', index=False)

if __name__ == "__main__":
    root_path = os.path.abspath(os.path.join(__file__, "../../.."))

    genotype = pd.read_csv(f'{root_path}/data/genotype_information/piQTL_genotype_matrix_dec2022.txt')
    PPI_list = pd.read_csv(f'{root_path}/data/pipeline/PPI_reference_barcodes.csv')['PPI']

    results_folder = f'{root_path}/results/'
    per_condition_folder = os.path.join(results_folder,'01_barcode_count/per_condition')


    for fld in [ os.path.join(results_folder, '01_barcode_count', 'per_PPI'), os.path.join(results_folder, '03_ppi_estimation'), os.path.join(results_folder, '03_ppi_estimation', 'logratio'), os.path.join(results_folder, '03_ppi_estimation', 'logratio', 'before_downsampling')]:
        try :
            os.mkdir(fld)
        except :
            print(f"{fld} already existing - Move to next steps")
    #### Split barcode count per PPI
    reads_per_ppi(per_condition_folder, PPI_list, results_folder)
    print("Estimate fitness for each PPI from barcode frequency")
    run_logratio_fitness_per_ppi(PPI_list, results_folder, root_path)
    print("Averaging of biological replicates ...")
    ppi_quantification_before_normalization(PPI_list, results_folder)
    print("Removing mean fitness of references strains ...")
    ppi_quantification_after_normalization(PPI_list, results_folder, genotype)
    print("Removing mean fitness of references strains & no PPI reference...")
    ppi_quantification_after_normalization_v2(PPI_list, results_folder, genotype)
