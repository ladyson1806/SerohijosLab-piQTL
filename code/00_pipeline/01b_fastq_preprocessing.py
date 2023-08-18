"""
NAME 01b_fastq_preprocessing.py

=========

DESCRIPTION

XXXXX

===========

INPUTS

XXXXX

============

OUTPUTS

XXXXX

=====

VERSION HISTORY

0.0.1   2022/12/11   Stable version.

=======

LICENCE

=======
2022, Copyright Savandara Besse (savandara.besse@umontreal.ca)


"""


from importlib.metadata import metadata
import os, re, subprocess
import numpy as np
import pandas as pd

from tqdm import tqdm

def get_env_cdts(x):
    if 'Time0' not in x :
        REP = x.split('_')[1]
        MTX = x.split('_')[2]
        DRUG = x.split('_')[3]
    else :
        REP = 'T0'
        MTX = 'noMTX'
        DRUG = 'noDrug'
    return f'{REP}_{MTX}_{DRUG}'

def demultiplexed_fastq(plate_folder, barcode_list):
    os.chdir(plate_folder)
    sample_metadata = plate_folder.split('/')[-1]
    print(f'Preprocessing reads for {sample_metadata}\n\n')

    sample_id_out = get_env_cdts(sample_metadata)
    adapter_removal_cmd = f'AdapterRemoval --file1 *_R1.fastq.gz --file2 *_R2.fastq.gz --basename {sample_id_out} --barcode-list {barcode_list} --barcode-mm-r1 2 --barcode-mm-r2 2  --threads 24 --gzip'
    subprocess.run(adapter_removal_cmd, shell=True)

    #### Stats on read lost
    init_read_count = int(subprocess.run('zcat *_R1.fastq.gz | grep "@" | wc -l', stdout=subprocess.PIPE,  shell=True).stdout.decode('utf-8').strip())
    lost_read_count = int(subprocess.run(f'zcat {sample_id_out}.unidentified_1.gz | grep "@" | wc -l', stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip())

    lost_read_ratio = (lost_read_count / init_read_count) * 100
    print(f'\nNumber of lost reads for {sample_id_out}: {lost_read_ratio} \n\n')

    subprocess.run(f'mv {sample_id_out}*.pair*.truncated.gz ../../demultiplexed/', shell=True)
    return sample_id_out


def merge_reads(sample_fastq):
    subprocess.run(f'pear -f {sample_fastq}.pair1.truncated.gz -r {sample_fastq}.pair2.truncated.gz -o ../merged/{sample_fastq} -j 24', shell=True)

if __name__ == "__main__":
    table_folder = ('../../piQTL/data/pipeline')
    PPI_list = pd.read_csv(os.path.join(table_folder, 'PPI_reference_barcodes.csv'))['PPI']
    barcode_list = os.path.join(table_folder,'PPI_reference_barcodes.txt')

    base_folder = '/media/UTELifeNAS/homes/Savvy/piQTL/DATA' #### Folder where are the final dataset with the concatanated fastq.gz are

    print('piQTL pipeline - Barcode extraction')

    ### 1. Demultiplexing + inner barcode trimming
    original_folder = os.path.join(base_folder,'original') #'/home/savvy/PROJECTS/PHD/DATA/original/'
    print('Step 1 - Demultiplexing + inner barcode trimming')
    for folder in sorted(os.listdir(original_folder)):
        for fold in ['demultiplexed', 'merged', 'pooled']:
            if not os.path.isdir(os.path.join(base_folder, fold)):
                os.mkdir(os.path.join(base_folder, fold))

        curr_folder = os.path.join(original_folder, folder)
        print(f'Woking on {folder} ...')
        print('Step 1 - Demultiplexing + inner barcode trimming')
        sample_name = demultiplexed_fastq(curr_folder, barcode_list)
        if not os.path.isdir(os.path.join(base_folder, sample_name)):
            os.mkdir(os.path.join(base_folder, sample_name))

        #### 2. Pair-end read merging files added in pooled folder
        print('Step 2 - Pair-end read merging')
        demultiplex_folder = os.path.join(base_folder,'demultiplexed') #'/home/savvy/PROJECTS/PHD/DATA/demultiplexed/'
        u_cdt_name = [ f'{sample_name}.{PPI}' for PPI in PPI_list ]

        print(u_cdt_name[0])
        os.chdir(demultiplex_folder)
        for condition in tqdm(np.unique(u_cdt_name)):
            merge_reads(condition)
        merged_folder = os.path.join(base_folder,'merged') #'/home/savvy/PROJECTS/PHD/DATA/merged/'
        os.chdir(merged_folder)
        subprocess.run(f'mv *.assembled.fastq ../pooled/', shell=True)

        ##### 3. Extract the barcodes
#         print('Step 3 - Barcode extraction')
#         pooled_folder = os.path.join(base_folder,'pooled') #'/home/savvy/PROJECTS/PHD/DATA/pooled/'
#         os.chdir(pooled_folder)
#         for condition in tqdm(np.unique(u_cdt_name)):
#             ##### With new barcode pattern, no need to have the cut step
#             bartender_extract_cmd = f'bartender_extractor_com -f {condition}.assembled.fastq -o {condition} -q "?" -p "TGGGC[5]CAGGTCTGAAGCTGTCGCAC[5]GAAAT" -m 4 -d both >> {condition}.log'
#             subprocess.run(bartender_extract_cmd, shell=True)

        #### 3. Move processed folders to associated condition folder and the final output
        os.chdir(base_folder)
        os.mkdir(os.path.join(base_folder,'barcode_extracted'))
        subprocess.run('rm -r demultiplexed merged', shell=True)
        subprocess.run(f'mv demultiplexed merged pooled {sample_name}', shell=True)
        subprocess.run(f'mv {sample_name} barcode_extracted')
