"""
NAME 01c_barcode_extraction.py

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

0.0.1   2022/12/12   Stable version.

=======

LICENCE

=======
2022, Copyright Savandara Besse (savandara.besse@umontreal.ca)


"""

import multiprocessing, os, signal, subprocess
from tqdm import tqdm

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)

def barcode_extraction(bartender_input):
    fastq = bartender_input
    output = f'{bartender_input.replace(".assembled.fastq","")}'
    log = f'{bartender_input.replace(".assembled.fastq","")}.log'
    bartender_extract_cmd = f'bartender_extractor_com -f {fastq} -o {output} -q "?" -p "TGGGC[5]CAGGTCTGAAGCTGTCGCAC[5]GAAAT" -m 4 -d both >> {log}'
    subprocess.run(bartender_extract_cmd, shell=True)


def main(arg):
    bartender_input = arg

    barcode_extraction(bartender_input)


if __name__ == "__main__":
    base_folder = '/home/savvy/PROJECTS/PHD/DATA/FINAL_DATASET'
    print('Step 4 - Barcode extraction')
    barcode_extracted_folder = os.path.join(base_folder,'barcode_extracted')
    for condition in os.listdir(barcode_extracted_folder):
        print(f"Barcode extraction for {condition}")
        pooled_folder = f'/home/savvy/PROJECTS/PHD/DATA/FINAL_DATASET/barcode_extracted/{condition}/pooled'
        args = [ os.path.join(pooled_folder, fastq) for fastq in tqdm(os.listdir(pooled_folder))]
        p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=10)
        try:
            pbar = tqdm(args, maxinterval=1.0, miniters=1, desc="Barcode extracted: ", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
            for _, result in enumerate(p.imap_unordered(main, args, chunksize=1)):
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
