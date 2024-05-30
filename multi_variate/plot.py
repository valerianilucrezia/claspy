import sys
sys.path.append('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy')

from queue import PriorityQueue
from claspy.utils import check_input_time_series, check_excl_radius
from claspy.clasp import ClaSPEnsemble
from claspy.validation import map_validation_tests
import numpy as np
import matplotlib.pyplot as plt

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_data(snv):
    sns.set_theme(style="white", font_scale=1.5)
    fig, axes = plt.subplots(2, 3, figsize=(25, 6))

    baf_plot = sns.scatterplot(data=snv, x="pos", y="median_baf", s=20, ax=axes[0,0], hue="cna_id", legend=False)
    baf_plot.axhline(0.5)
    baf_plot_line = sns.lineplot(data=snv, x = 'pos', y="median_baf", ax=axes[1,0], legend=False)
    

    dr_plot = sns.scatterplot(data=snv, x="pos", y="median_dr", s=20, ax=axes[0,1], hue="cna_id",  legend=False)
    dr_plot.axhline(1)
    dr_plot_line = sns.lineplot(data=snv, x = 'pos', y="median_dr", ax=axes[1,1], legend=False)

    vaf_plot = sns.scatterplot(data=snv, x="pos", y="vaf", s=20, ax=axes[0,2], hue="cna_id",  legend=False)
    vaf_plot_line = sns.lineplot(data=snv, x = 'pos', y="vaf", ax=axes[1,2], legend=False)

    axes[0,0].set_ylim(0,1)
    axes[0,2].set_ylim(0,1)
    return


def get_original_data_plot(data = '/Users/lucreziavaleriani/Desktop/seq_res_snv_smooth.csv', 
                           save = False, 
                           out_file = '/Users/lucreziavaleriani/Desktop/original.png'):
    
    snv = pd.read_csv(data)
    snv.fillna(value='normal', axis = 0, inplace = True)
    snv.sort_values(by = ['pos'], inplace = True)
    snv['id'] = [i for i in range(snv.shape[0])]
                    
    bps_max = list(snv.groupby(['cna_id']).max(['pos']).id)
    bps_min = list(snv.groupby(['cna_id']).min(['pos']).id)
    bps = np.array(bps_max + bps_min)                
    vaf_bps = snv.groupby('cna_id')
                    
    vaf = np.array(snv.vaf)
    baf = np.array(snv.median_baf)
    dr = np.array(snv.median_dr)
                
    plot_data(snv)            
    if save and out_file!= '':
        plt.savefig(fname = os.path.join(out_file, 'original_data.png'))
    plt.close()
    
    return vaf, baf, dr, bps


def plot_profile(dr_profile, baf_profile, vaf_profile, 
                 true_bps, pred_bps = [],
                 title = '', 
                 save = False, 
                 out_file = '', 
                 mode = ''):
    
    # fig = plt.figure(figsize=(10, 6))
    fig, axs = plt.subplots(nrows = 3, ncols = 1, figsize=(10, 8))
    # ax = plt.subplot(111)

    axs[0].plot(dr_profile, label = "DR")
    axs[1].plot(baf_profile, label = 'BAF')
    axs[2].plot(vaf_profile, label = 'VAF')
    
    axs[0].set_ylabel('DR')
    axs[1].set_ylabel('BAF')
    axs[2].set_ylabel('VAF')
    
    axs[0].vlines(true_bps, ymin = min(dr_profile) - 0.05, ymax = max(dr_profile) + 0.05, colors = 'tab:green', label = 'True BP', linestyles = 'dashed')
    axs[0].vlines(pred_bps, ymin = min(dr_profile), ymax = max(dr_profile), colors = 'tab:olive',  label = 'Predicted BP')

    axs[1].vlines(true_bps, ymin = min(baf_profile) - 0.05, ymax = max(baf_profile) + 0.05, colors = 'tab:green', label = 'True BP', linestyles = 'dashed')
    axs[1].vlines(pred_bps, ymin = min(baf_profile), ymax = max(baf_profile), colors = 'tab:olive', label = 'Predicted BP')
    
    axs[2].vlines(true_bps, ymin = min(vaf_profile) - 0.05, ymax = max(vaf_profile) + 0.05, colors = 'tab:green', label = 'True BP', linestyles = 'dashed')
    axs[2].vlines(pred_bps, ymin = min(vaf_profile), ymax = max(vaf_profile), colors = 'tab:olive', label = 'Predicted BP')


    #axs.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
    fig.suptitle(title)
    fig.tight_layout()
    
    if save and out_file != '':
        plt.savefig(os.path.join(out_file, mode + '_result.png'), dpi = 500)
    plt.close()
    return