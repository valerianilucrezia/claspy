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
    
    return vaf, baf, dr, bps


def plot_profile(dr_profile, baf_profile, vaf_profile, 
                 true_bps, pred_bps = [],
                 title = '', 
                 save = False, 
                 out_file = '', 
                 mode = ''):
    
    fig = plt.figure(figsize=(10, 6))
    ax = plt.subplot(111)

    plt.plot(dr_profile, label = "DR")
    plt.plot(baf_profile, label = 'BAF')
    plt.plot(vaf_profile, label = 'VAF')
    
    plt.vlines(true_bps, ymin = 0, ymax = 2, colors = 'tab:pink', label = 'True BP')
    plt.vlines(pred_bps, ymin = 0, ymax = 2, colors = 'tab:cyan', label = 'Predicted BP')

    ax.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
    plt.title(title)
    plt.tight_layout()
    #plt.show()
    
    if save and out_file != '':
        plt.savefig(os.path.join(out_file, mode + '_result.png'))
    
    return