import cv2, os
import numpy as np 
import pandas as pd
import seaborn as sns 
import matplotlib as mpl
import matplotlib.pyplot as plt

from tqdm import tqdm 
from functools import reduce
from assocplots.manhattan import *
from assocplots.qqplot  import *
from matplotlib.lines import Line2D


# set matplotlib default parameters
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = "Helvetica"
mpl.rcParams["font.size"] = 4
mpl.rcParams["axes.titlesize"] = 4
mpl.rcParams["xtick.labelsize"] = 4
mpl.rcParams["xtick.major.size"] = 2
mpl.rcParams["xtick.minor.size"] = 1.2
mpl.rcParams["ytick.labelsize"] = 4
mpl.rcParams["ytick.major.size"] = 2
mpl.rcParams["ytick.minor.size"] = 1.2
mpl.rcParams["lines.linewidth"] = 1.2
mpl.rcParams["lines.markersize"] = 2
mpl.rcParams["lines.markeredgewidth"] = 0.5
mpl.rcParams["boxplot.flierprops.markersize"] = 2
mpl.rcParams["boxplot.flierprops.markeredgewidth"] = 0.5
mpl.rcParams["legend.fontsize"] = 4

CM = 1/2.54 # centimeters in inches
QUARTER = 4.5 * CM
THIRD = 6 * CM
HALF = 9 * CM
WHOLE = 18 * CM


def panel_1DV(PPI, DRUG):
    i = 0
    for MTX in ['MTX', 'noMTX']:
        REPS = []
        LT = []
        for REP in ['A', 'B']:
            PPI_table = pd.read_csv(f'../results/01_barcode_count/per_PPI/{PPI}_read_number.csv')
            R = PPI_table[['strain_id', f'T0_noMTX_noDrug.{PPI}', f'{REP}_{MTX}_{DRUG}.{PPI}']]
            R[f'{REP}_{MTX}_{DRUG}.{PPI}'] = R[f'{REP}_{MTX}_{DRUG}.{PPI}'] + 1 
            R[f'T0.{PPI}'] = R[f'T0_noMTX_noDrug.{PPI}'] + 1 

            R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] = (R[f'{REP}_{MTX}_{DRUG}.{PPI}'] / sum(R[f'{REP}_{MTX}_{DRUG}.{PPI}'])) * 10**6
            R[f'T0.{PPI}_RPM'] = (R[f'T0.{PPI}'] / sum(R[f'T0.{PPI}'])) * 10**6
            R['logratio_Fitness'] = np.log2(R[f'{REP}_{MTX}_{DRUG}.{PPI}_RPM'] / R[f'T0.{PPI}_RPM'])
            RPM_cols = [ col for col in R.columns if 'RPM' in col ]
            RPM_cols.insert(0, 'strain_id')
            RPM_cols.insert(-1, 'logratio_Fitness')
            #  display(R[RPM_cols])
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
        ref_colors = ['red', 'maroon', 'cyan', 'darkblue'] # 'paleturquoise', 'darkturquoise', 'darkcyan', 'yellow', 'orange', 'brown'
 
        g, g_axes = plt.subplots(ncols=2, nrows=2, figsize=(QUARTER, QUARTER), sharex=False, sharey=False)
        custom_lines = [
            Line2D([0], [0], color='red', linestyle='--', lw=2),
            Line2D([0], [0], color='maroon', linestyle='--', lw=2),
            Line2D([0], [0], color='cyan', linestyle='--', lw=2),
            Line2D([0], [0], color='marine', linestyle='--', lw=2),
            # Line2D([0], [0], color='orange', linestyle='--', lw=2),
            # Line2D([0], [0], color='brown', linestyle='--', lw=2),
            ]
        
        # g_axes[0].legend(custom_lines, ['43','43_ref1', '43_ref2', '599', '599_ref1', '599_ref2'], loc='center right')
        # g_axes[1].legend(custom_lines, ['43','43_ref1', '43_ref2', '599', '599_ref1', '599_ref2'], loc='center right', fontsize=8)


        A_segregants = LT[0][~LT[0]['strain_id'].str.contains('ref')].reset_index(drop=True)
        A_ref_segregants = LT[0][LT[0]['strain_id'].str.contains('ref') ].reset_index(drop=True) # | (scatter_table['strain_id'].isin(['43','599']))
        

        sns.lineplot(data=A_segregants, x='Time', y='RPM', hue='logratio_Fitness', ax=g_axes[0,0], legend=False)
        sns.lineplot(data=A_ref_segregants, x='Time', y='RPM', hue='strain_id',  palette=ref_colors, linestyle='dashed', ax=g_axes[0,0], legend=False)
        g_axes[0,0].set_yscale('log')
        xticks = [0,96]
        xticklabels = ['0', '96']
        g_axes[0,0].set_xticks(xticks)
        g_axes[0,0].set_xticklabels(xticklabels)
        g_axes[0,0].set_xlabel('Time (H)')
        # g_axes[0,0].set_title('Replicate A')

        sns.kdeplot(data=LT[0]['logratio_Fitness'], color='black', ax=g_axes[0,1])
        # g_axes[0,1].set_title('Replicate A')
        g_axes[0,1].set_xlabel('Fitness')
        
        k1 = 0
        for ref in np.unique(A_ref_segregants['strain_id']):
            ref_fitness = A_ref_segregants[A_ref_segregants['strain_id'] == ref]['logratio_Fitness'].values[0]
            # print(ref, ref_fitness)
            g_axes[0,1].axvline(x=ref_fitness, linestyle='dashed', color=ref_colors[k1])
            k1 += 1

        plt.tight_layout()

        B_segregants = LT[1][~LT[1]['strain_id'].str.contains('ref')].reset_index(drop=True)
        B_ref_segregants = LT[1][LT[1]['strain_id'].str.contains('ref') ].reset_index(drop=True)
        
        sns.lineplot(data=B_segregants, x='Time', y='RPM', hue='logratio_Fitness', legend=False, ax=g_axes[1, 0])
        sns.lineplot(data=B_ref_segregants, x='Time', y='RPM', hue='strain_id', palette=ref_colors, linestyle='dashed', legend=False, ax=g_axes[1, 0])
        g_axes[1,0].set_yscale('log')
        g_axes[1,0].set_xticks(xticks)
        g_axes[1,0].set_xticklabels(xticklabels)
        g_axes[1,0].set_xlabel('Time (H)')
        # g_axes[1,0].set_title('Replicate B')

        sns.kdeplot(data=LT[1]['logratio_Fitness'], color='black', ax=g_axes[1, 1])
        # g_axes[1,1].set_title('Replicate B')
        g_axes[1,1].set_xlabel('Fitness')

        k2 = 0
        for ref in np.unique(B_ref_segregants['strain_id']):
            ref_fitness = B_ref_segregants[B_ref_segregants['strain_id'] == ref]['logratio_Fitness'].values[0]
            # print(ref, ref_fitness)
            g_axes[1,1].axvline(x=ref_fitness, linestyle='dashed', color=ref_colors[k2])
            k2 += 1
        
        # plt.suptitle(f'{PPI} under {DRUG} ({MTX})')

        g.savefig('/home/savvy/PROJECTS/PHD/piQTL/manuscript/figures/FIGURE_1/FIGURE_1dv.eps', dpi=300)

def manhattanplot_graph(QTL, label):
    f = plt.figure(figsize=(QUARTER,QUARTER/3))
    manhattan(QTL['P'], QTL['BP'], QTL['CHR'], '',
    plot_type='single',
    chrs_plot=None,
    chrs_names=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', 'M'],
    cut = 0,
    title=label,
    xlabel='Chromosome',
    ylabel='-log10(p-value)',
    lines=[3],
    lines_colors=['black'],
    lines_styles=['--'],
    scaling = '-log10') 

    f.savefig('/home/savvy/PROJECTS/PHD/piQTL/manuscript/figures/FIGURE_1/FIGURE_1dvia.eps', dpi=300)

def qqplot_graph(QTL, label) :
    f = plt.figure(figsize=(QUARTER/2,QUARTER/2)) 
    ax = plt.subplot(111)
    qqplot([QTL['P']], 
        [label], 
        color=['grey'], 
        alpha=0.95,
        fill_dens=[0.2], 
        error_type='theoretical', 
        distribution='beta',
        title='')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5))

    f.savefig('/home/savvy/PROJECTS/PHD/piQTL/manuscript/figures/FIGURE_1/FIGURE_1dvib.eps', dpi=300)

if __name__ == "__main__":
    PPI_list = pd.read_csv('../data/pipeline/PPI_reference_barcodes.csv')['PPI']  
    results_folder = '/home/savvy/PROJECTS/PHD/piQTL/results'
    
    PPI = 'ERG11_PIS1'
    DRUG = 'Fluconazole'
    
    # panel_1DV(PPI, DRUG)


    table_path = '/home/savvy/PROJECTS/PHD/piQTL/data/QTL'
    rMVP_qtl_res = pd.read_csv(os.path.join(table_path, f'{PPI}_MTX_{DRUG}_avg_logratio_Fitness_minus_ref.csv'))
    print(rMVP_qtl_res.columns)
    manhattanplot_graph(rMVP_qtl_res, '')
    qqplot_graph(rMVP_qtl_res, '')