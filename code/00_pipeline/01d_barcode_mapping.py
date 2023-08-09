"""
NAME 01d_barcode_mapping.py

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
import multiprocessing, os, signal, textdistance
import numpy as np
import pandas as pd 
import seaborn as sns 
import scipy.stats as stats 
import matplotlib.pyplot as plt  

from tqdm import tqdm
from functools import reduce
tqdm.pandas()

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)


def retrieve_strain_id(x, barcode_reference_library):
    all_true_barcodes = barcode_reference_library['revComp'].values
    distances = []
    for true_barcode in all_true_barcodes :
        distances.append(textdistance.hamming.distance(x,true_barcode))
    distances = np.array(distances)
    closest_barcode = sorted([ (distances[i], all_true_barcodes[i]) for i in range(len(distances)) if distances[i] <= 2], reverse=True)
    if len(closest_barcode) == 1:
        candidate_barcode = closest_barcode[0][1]
        return barcode_reference_library[barcode_reference_library['revComp'] == candidate_barcode]['strain_number'].values[0]
    if len(closest_barcode) > 1 :
        if closest_barcode[-1][0] == 1 :
            candidate_barcode = closest_barcode[-1][1]
            return barcode_reference_library[barcode_reference_library['revComp'] == candidate_barcode]['strain_number'].values[0]
    return np.nan

def retrieve_strain_id_strict(x, all_true_barcodes):
    if x in all_true_barcodes.keys() : 
        return all_true_barcodes[x]
    else :
        for true_barcode in all_true_barcodes.keys():
            if textdistance.hamming.distance(x,true_barcode) == 1 :
                return all_true_barcodes[true_barcode]
        return np.nan

def barcode_mapping(CDT):
    ALL = []
    read_mapping_stat = []
    pooled_folder = f'/home/savvy/PROJECTS/PHD/DATA/FINAL_DATASET/barcode_extracted/{CDT}/pooled'
    for table in os.listdir(pooled_folder):
        if 'barcode.txt' in table :
            PPI = table.split('.')[-1].replace("_barcode.txt","")
            TMP = pd.read_csv(f'{pooled_folder}/{table}',names=['barcode','count'])
            total_reads = len(TMP)
            TMP['strain_id'] = TMP['barcode'].apply(retrieve_strain_id_strict, args=(all_true_barcodes,))
            # TMP['strain_id'] = TMP['barcode'].apply(retrieve_strain_id, args=(barcode_reference_library,))
            mapped_reads = sum(TMP.dropna().groupby('strain_id').count()['count'])
            print(table, total_reads, mapped_reads, len(TMP.dropna().groupby('strain_id').count()))
            retrieve_barcode = mapped_reads / total_reads * 100
            TMP = TMP.dropna().groupby('strain_id').count().reset_index().rename(columns={'count':f'{table.replace("_barcode.txt","")}'}).drop(columns=['barcode'])
            ALL.append(TMP)
            read_mapping_stat.append([f'{table.replace("_barcode.txt","")}', total_reads, mapped_reads, retrieve_barcode, len(TMP)])
    FINAL_TIME_MERGED = reduce(lambda  left,right: pd.merge(left,right,on=['strain_id'], how='outer'), ALL).replace(np.nan, 0)

    barcode_mapping_stats = pd.DataFrame(read_mapping_stat, columns=['Condition', 'total_reads', 'mapped_reads', '%_barcode_mapped', 'nb_strains'])
    barcode_mapping_stats.to_csv(f'../../results/01_barcode_count/log/{CDT}.log', index=False)
    FINAL_TIME_MERGED.to_csv(f'../../results/01_barcode_count/per_condition/{CDT}_read_number.csv', index=False)


def main(folder):
    barcode_mapping(folder)


if __name__ == "__main__":
    barcode_reference_library = pd.read_csv('../../data/pipeline/barcode_collection_with_ubarcodes.csv')
    barcode_reference_library = barcode_reference_library[['strain_number','revComp']]
    print(barcode_reference_library)
    all_true_barcodes = dict(zip(barcode_reference_library.revComp, barcode_reference_library.strain_number)) 
    
    barcode_extraction_folder = '/home/savvy/PROJECTS/PHD/DATA/FINAL_DATASET/barcode_extracted/'

    base_folder = '../../results'
    try :
        for fld in [os.path.join(base_folder, '01_barcode_count'), os.path.join(base_folder, '01_barcode_count', 'per_condition'), os.path.join(base_folder, '01_barcode_count', 'log')]:
            os.mkdir(fld)
    except :
        print("Folder already existing - Move to next steps")
    print('Step 5 - Barcode Mapping')
    myfolders =  [ folder for folder in sorted(os.listdir(barcode_extraction_folder)) ]
    p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=12)
    try:
        pbar = tqdm(myfolders, maxinterval=1.0, miniters=1, desc="Barcode mapped: ", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
        for _, result in enumerate(p.imap_unordered(main, myfolders, chunksize=1)):
            pbar.update(1)  # Everytime the iteration finishes, update the global progress bar

        pbar.close()
        p.close()
        p.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt, terminating workers.")
        pbar.close()
        p.terminate()
        p.join()
        exit(1)