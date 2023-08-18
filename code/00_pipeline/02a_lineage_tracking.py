import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import cv2, multiprocessing, os, signal
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

from tqdm import tqdm
from functools import reduce
from IPython.display import display
from matplotlib.lines import Line2D


plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams.update({'font.size': 7})

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)

def vconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    # take minimum width
    w_min = min(img.shape[1] for img in img_list)

    # resizing images
    im_list_resize = [cv2.resize(img, (w_min, int(img.shape[0] * w_min / img.shape[1])), interpolation = interpolation) for img in img_list]
    # return final image
    return cv2.vconcat(im_list_resize)

def hconcat_resize(img_list, interpolation=cv2.INTER_CUBIC):
    # take minimum hights
    h_min = min(img.shape[0] for img in img_list)

    # image resizing
    im_list_resize = [cv2.resize(img, (int(img.shape[1] * h_min / img.shape[0]), h_min), interpolation = interpolation)
                      for img in img_list]

    # return final image
    return cv2.hconcat(im_list_resize)

def fitness_comparison(PPI, DRUG, res, root):
    f, axes = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(12,8))
    i = 0
    REPS = []
    LT = []
    for REP in ['A', 'B']:
        for MTX in ['MTX', 'noMTX']:
            PPI_table = pd.read_csv(f'{root}/results/01_barcode_count/per_PPI/{PPI}_read_number.csv')
            R = PPI_table[['strain_id', f'T0_noMTX_noDrug.{PPI}', f'{REP}_{MTX}_{DRUG}.{PPI}']]
            R[f'{REP}_{MTX}_{DRUG}.{PPI}'] = R[f'{REP}_{MTX}_{DRUG}.{PPI}'] + 1
            R[f'T0.{PPI}'] = R[f'T0_noMTX_noDrug.{PPI}'] + 1

            R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] = (R[f'{REP}_{MTX}_{DRUG}.{PPI}'] / sum(R[f'{REP}_{MTX}_{DRUG}.{PPI}'])) * 10**6
            R[f'T0.{PPI}_RPM'] = (R[f'T0.{PPI}'] / sum(R[f'T0.{PPI}'])) * 10**6
            R['logratio_Fitness'] = np.log2(R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] / R[f'T0.{PPI}_RPM'])
            RPM_cols = [ col for col in R.columns if 'RPM' in col ]
            RPM_cols.insert(0, 'strain_id')
            RPM_cols.insert(-1, 'logratio_Fitness')
            REPS.append(R[RPM_cols])

            table_per_timepoint = []
            TMP = R[['strain_id', f'T0.{PPI}_RPM', f'logratio_Fitness']].rename(columns={f'T0.{PPI}_RPM': 'RPM'})
            TMP['Time'] = '0'
            table_per_timepoint.append(TMP)
            TMP = R[['strain_id', f'{REP}_{MTX}_{DRUG}.{PPI}_RPM', f'logratio_Fitness']].rename(columns={f'{REP}_{MTX}_{DRUG}.{PPI}_RPM': 'RPM'})
            TMP['Time'] = '96'
            table_per_timepoint.append(TMP)

            scatter_table =  pd.concat(table_per_timepoint).reset_index(drop=True)
            scatter_table['Time'] = scatter_table['Time'].astype(int)
            LT.append(scatter_table)

        #### Lineage tracking plot
        ref_colors = ['paleturquoise', 'darkturquoise', 'darkcyan', 'yellow', 'orange', 'brown']

        g, g_axes = plt.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(12,12))

        custom_lines = [
            Line2D([0], [0], color='paleturquoise', linestyle='--', lw=2),
            Line2D([0], [0], color='darkcyan', linestyle='--', lw=2),
            Line2D([0], [0], color='darkturquoise', linestyle='--', lw=2),
            Line2D([0], [0], color='yellow', linestyle='--', lw=2),
            Line2D([0], [0], color='orange', linestyle='--', lw=2),
            Line2D([0], [0], color='brown', linestyle='--', lw=2),
            ]

        g_axes[0,1].legend(custom_lines, ['43','43_ref1', '43_ref2', '599', '599_ref1', '599_ref2'], loc='center right', fontsize=8)
        g_axes[1,1].legend(custom_lines, ['43','43_ref1', '43_ref2', '599', '599_ref1', '599_ref2'], loc='center right', fontsize=8)


        A_segregants = LT[0][~LT[0]['strain_id'].str.contains('ref')].reset_index(drop=True)
        A_ref_segregants = LT[0][LT[0]['strain_id'].str.contains('ref') | (scatter_table['strain_id'].isin(['43','599'])) ].reset_index(drop=True)

        sns.lineplot(data=A_segregants, x='Time', y='RPM', hue='logratio_Fitness', legend=False, ax=g_axes[0, 0])
        sns.lineplot(data=A_ref_segregants, x='Time', y='RPM', hue='strain_id',  palette=ref_colors, linestyle='dashed', legend=False, ax=g_axes[0, 0])
        g_axes[0, 0].set_yscale('log')
        g_axes[0, 0].set_title(f'{PPI}_{DRUG}_{MTX} - Replicate A')
        sns.kdeplot(data=LT[0]['logratio_Fitness'], ax=g_axes[0, 1])
        g_axes[0, 1].set_title(f'{PPI}_{DRUG}_{MTX} - Replicate A')

        k1 = 0
        for ref in np.unique(A_ref_segregants['strain_id']):
            ref_fitness = A_ref_segregants[A_ref_segregants['strain_id'] == ref]['logratio_Fitness'].values[0]
            # print(ref, ref_fitness)
            g_axes[0,1].axvline(x=ref_fitness, linestyle='dashed', color=ref_colors[k1])
            k1 += 1

        B_segregants = LT[1][~LT[1]['strain_id'].str.contains('ref')].reset_index(drop=True)
        B_ref_segregants = LT[1][LT[1]['strain_id'].str.contains('ref') | (scatter_table['strain_id'].isin(['43','599'])) ].reset_index(drop=True)

        sns.lineplot(data=B_segregants, x='Time', y='RPM', hue='logratio_Fitness', legend=False, ax=g_axes[1, 0])
        sns.lineplot(data=B_ref_segregants, x='Time', y='RPM', hue='strain_id', palette=ref_colors, linestyle='dashed', legend=False, ax=g_axes[1, 0])
        g_axes[1, 0].set_yscale('log')
        g_axes[1, 0].set_title(f'{PPI}_{DRUG}_{MTX} - Replicate B')
        sns.kdeplot(data=LT[1]['logratio_Fitness'], ax=g_axes[1, 1])
        g_axes[1, 1].set_title(f'{PPI}_{DRUG}_{MTX} - Replicate B')

        k2 = 0
        for ref in np.unique(B_ref_segregants['strain_id']):
            ref_fitness = B_ref_segregants[B_ref_segregants['strain_id'] == ref]['logratio_Fitness'].values[0]
            # print(ref, ref_fitness)
            g_axes[1,1].axvline(x=ref_fitness, linestyle='dashed', color=ref_colors[k2])
            k2 += 1

        g.savefig(f'./tmp_{MTX}_bd.png', dpi=300)

    #### Fitness comparison
    REP_MERGED = pd.merge(REPS[0], REPS[1], on=['strain_id', f'T0.{PPI}_RPM'], suffixes=['_A', '_B'])
    REP_MERGED_REF = REP_MERGED[REP_MERGED['strain_id'].str.contains('_ref')]
    # print(REP_MERGED)
    REP_MERGED['logratio_Fitness_A'] = REP_MERGED['logratio_Fitness_A'] - REP_MERGED_REF['logratio_Fitness_A'].mean()
    REP_MERGED['logratio_Fitness_B'] = REP_MERGED['logratio_Fitness_B'] - REP_MERGED_REF['logratio_Fitness_B'].mean()
    REP_MERGED['MTX_condition'] = f'{REP}_{PPI}_{DRUG}_{MTX}'
    # print(REP_MERGED)
    # print(REP_MERGED)
    # REP_MERGED = REP_MERGED[REP_MERGED['strain_id'].isin(['17', '17_2', '40', '42_2', '180', '180_2', '599_ref1', '599_ref2', '43_ref1', '43_ref2'])]
    sns.scatterplot(data=REP_MERGED, x='logratio_Fitness_A', y='logratio_Fitness_B', color='black', ax=axes[i])
    axes[i].set_xlabel('Fitness (Replicate A)')
    axes[i].set_ylabel('Fitness (Replicate B)')
    sns.scatterplot(data=REP_MERGED[REP_MERGED['strain_id'].str.contains('_ref')], x='logratio_Fitness_A', y='logratio_Fitness_B', color='red', ax=axes[i])
    # axes[i].axline([-1,-1],[1,1], linestyle='dashed', color='grey')
    corr, pval = stats.pearsonr(REP_MERGED['logratio_Fitness_A'], REP_MERGED['logratio_Fitness_B'])
    axes[i].annotate(f"spearman's corr: {round(corr,3)}\npval: {'{:.2e}'.format(pval)}", xy=[0, -60], xycoords="axes points")
    axes[i].set_title(f'{PPI} under {DRUG} ({MTX})'.replace('_',':'))
    i += 1

    res[f'{PPI}_{DRUG}_{MTX}'] = {}
    res[f'{PPI}_{DRUG}_{MTX}']['corr'] = corr
    res[f'{PPI}_{DRUG}_{MTX}']['pval'] = pval
    f.savefig('./tmp_fitness_bd.png', dpi=300)
    plt.close()

    img1 = cv2.imread('./tmp_noMTX_bd.png')
    img2 = cv2.imread('./tmp_MTX_bd.png')
    img3 = cv2.imread('./tmp_fitness_bd.png')

    img_v_resize = vconcat_resize([img1, img2, img3])
    cv2.imwrite(f'{root}/results/02_lineage_tracking/{PPI}_{DRUG}.png', img_v_resize)
    return REP_MERGED[['strain_id', 'logratio_Fitness_A', 'logratio_Fitness_B', 'MTX_condition']], res, f'{root}/results/02_lineage_tracking/{PPI}_{DRUG}.png'

def merge_plots(PPI, ALL_DRUGS):
    all_imgs = [ cv2.imread(img) for img in ALL_DRUGS ]
    img_h_resize = hconcat_resize(all_imgs)
    cv2.imwrite(f'../../results/02_lineage_tracking/merged/{PPI}.png', img_h_resize)
    remove_imgs = [ os.remove(img) for img in ALL_DRUGS ]

def main(args):
    PPI = args[0]
    ALL_DRUGS = args[1]
    merge_plots(args)

if __name__ == "__main__":
    root_path = os.path.abspath(os.path.join(__file__, "../../.."))

    PPI_list = pd.read_csv(f'{root_path}/data/pipeline/PPI_reference_barcodes.csv')['PPI']
    results_folder = f'{root_path}/results'
    try :
        for fld in [os.path.join(results_folder, '02_lineage_tracking'),  os.path.join(results_folder, '02_lineage_tracking', 'before_downsampling')]:
            os.mkdir(fld)
    except :
        print("Folder already existing - Move to next steps")
    args = []
    ALL_TABLE = []
    res = {}
    # for PPI in tqdm(PPI_list) :
    #     ALL_DRUGS = []
    #     print(f'Lineage Tracking for {PPI} in progress')
    #     for MTX in ['noMTX', 'MTX']:
    #         for DRUG in ['noDrug', '5-FC', 'Fluconazole', 'Metformin', 'Trifluoperazine'] :
    #             TABLE, res = fitness_comparison(PPI, DRUG, MTX, res)
    #             ALL_TABLE.append(TABLE)

    #### Dump all fitness of A and B in a master table
    # pd.concat(ALL_TABLE).to_csv('../01_experimental_checks/A_vs_B_annotated_replicates.csv', index=False)

    PPI_list_10 = ['ALO1_ADE17', 'ERG11_MID2'] # 'ERG11_MID2', 'SRO9_SGN1', 'ERG11_ORM2', 'HNM1_MID2', 'ERV25_ORM2', 'ERG11_TPO4', 'ERG11_PIS1', 'GTR1_SLM4', 'HNM1_ERG11'
    for PPI in tqdm(PPI_list_10) :
        ALL_DRUGS = []
        print(f'Lineage Tracking for {PPI} in progress')
        for DRUG in ['noDrug', '5-FC', 'Fluconazole', 'Metformin', 'Trifluoperazine'] :
            TABLE, res , drug_image = fitness_comparison(PPI, DRUG, res, root_path)
            ALL_DRUGS.append(drug_image)
            ALL_TABLE.append(TABLE)

    # #### Dump all fitness of A and B in a master table
    # pd.concat(ALL_TABLE).to_csv('./A_vs_B_10.csv', index=False)

    # f = open('./lineage_tracking_stats.csv', 'w')
    # f.write('Condition,Correlation,p-value\n')
    # for PPI in tqdm(res.keys()) :
    #     f.write(f"{PPI},{res[PPI]['corr']},{res[PPI]['pval']}\n")
    # f.close()


    #### Merge plots
    args.append((PPI, ALL_DRUGS))
    p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=12)
    try:
        pbar = tqdm(args, maxinterval=1.0, miniters=1, desc="Barcode mapped: ", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
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
