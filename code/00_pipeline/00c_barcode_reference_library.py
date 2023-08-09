#!/usr/bin/env python3

"""
NAME 00_barcode_reference_library.py

=========

DESCRIPTION 

Build the barcode table (Well Index) needed by Adapter Removal based on the plate design 

===========

INPUTS

- CSV file : ../../data/pipeline/sample_position_and_barcodes.csv

============

OUTPUTS

- CSV file : ../../data/pipeline/PPI_reference_barcodes.csv #### Modified version of sample_position_and_barcodes.csv 
                                                                with addition of 3 columns : N_barcode, S_barcode, PPI
- CSV file : ../../data/pipeline/PPI_reference_barcodes.txt #### Required file needed for 01_barcode_extraction.py


=====

VERSION HISTORY

0.0.1   2022/12/11   Stable version.

=======

LICENCE

=======
2022, Copyright Savandara Besse (savandara.besse@umontreal.ca)

"""

import pandas as pd 



#### Nseries - Left inner barcode (already reverse complement)
N_series = {
    'N701': 'TCGCCTTA', 
    'N702':	'CTAGTACG',
    'N703':	'TTCTGCCT',
    'N704':	'GCTCAGGA',
    'N705':	'AGGAGTCC',
    'N706':	'CATGCCTA',
    'N707':	'GTAGAGAG',
    'N708':	'CCTCTCTG',
    'N709':	'AGCGTAGC',
    'N710':	'CAGCCTCG',
    'N711':	'TGCCTCTT',
    'N712':	'TCCTCTAC'
}

#### Sseries - Right inner barcode (already reverse complement)
S_series = {
    'S501':	'TAGATCGC',
    'S502':	'CTCTCTAT',
    'S503':	'TATCCTCT',
    'S504':	'AGAGTAGA',
    'S505':	'GTAAGGAG',
    'S506':	'ACTGCATA',
    'S507':	'AAGGAGTA',
    'S508':	'CTAAGCCT',
    'S510':	'CGTCTAAT',
    'S511':	'TCTCTCCG',
    'S513':	'TCGACTAG',
    'S515':	'TTCTAGCT'

}

def get_reference_barcode(x, series_type):
    if 'N' in series_type :
        return N_series[x]
    if 'S' in series_type :
        return S_series[x]

if __name__ == "__main__":
    plate_design = pd.read_csv('../../data/pipeline/sample_position_and_barcodes.csv')
    plate_design['N_barcode'] = plate_design['N-index'].apply(get_reference_barcode, args=('N',))
    plate_design['S_barcode'] = plate_design['S-index'].apply(get_reference_barcode, args=('S',))
    plate_design['PPI'] = [ plate_design['ppi'][i].replace(':','_') for i in plate_design.index ]
    print(plate_design)

    plate_design.to_csv('../../data/pipeline/PPI_reference_barcodes.csv', index=False)
    plate_design[['PPI', 'N_barcode', 'S_barcode']].to_csv('../../data/pipeline/PPI_reference_barcodes.txt', index=False, sep='\t', header=False)
