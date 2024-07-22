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

# TODO: this is not called from what I can tell and is duplicated in get_original_data_plot
def get_data(data = '/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_nanopore/data/smooth_data.tsv'):
    snv = pd.read_csv(data, sep = '\t')
    snv.fillna(value=0, axis = 0, inplace = True)
    
    vaf = np.array(snv.gt_AF)
    baf = np.array(snv.median_baf)
    dr = np.array(snv.median_dr)
    maf = np.array(snv.median_meth)
    
    return vaf, baf, dr, maf


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


def get_original_data_plot(data = '', 
                           save = False, 
                           out_file = '/Users/lucreziavaleriani/Desktop/original.png'):
    # return dict to keep track of variable names
    ret_dict = {}
    snv = pd.read_csv(data)
    snv.fillna(value='normal', axis = 0, inplace = True)
    snv.sort_values(by = ['pos'], inplace = True)
    snv['id'] = [i for i in range(snv.shape[0])]
                    
    bps_max = list(snv.groupby(['cna_id']).max(['pos']).id)
    bps_min = list(snv.groupby(['cna_id']).min(['pos']).id)
    ret_dict["bps"] = np.array(bps_max + bps_min)                
    vaf_bps = snv.groupby('cna_id')
                    
    ret_dict["vaf"] = np.array(snv.vaf)
    ret_dict["baf"] = np.array(snv.median_baf)
    ret_dict["dr"] = np.array(snv.median_dr)
                
    plot_data(snv)            
    if save and out_file!= '':
        plt.savefig(fname = os.path.join(out_file, 'original_data.png'))
    plt.close()
    
    return ret_dict


def plot_profile(multivariate_clasp_objects, 
                 true_bps, pred_bps = [],
                 title = '', 
                 save = False, 
                 out_file = '', 
                 mode = ''):
    
    # fig = plt.figure(figsize=(10, 6))
    fig, axs = plt.subplots(nrows = 3, ncols = 1, figsize=(10, 8))
    # ax = plt.subplot(111)
    variables = multivariate_clasp_objects.keys()
    for i in range(0, len(variables)-1):
        profile = multivariate_clasp_objects[variables[i]].profile
        axs[i].plot(profile, label = variables)
        axs[i].set_ylabel(variables[i])
        axs[i].vlines(true_bps, ymin = min(profile) - 0.05, ymax = max(profile) + 0.05, colors = 'tab:green', label = 'True BP', linestyles = 'dashed')
        axs[i].vlines(pred_bps, ymin = min(profile), ymax = max(profile), colors = 'tab:olive',  label = 'Predicted BP')


    #axs.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
    fig.suptitle(title)
    fig.tight_layout()
    
    if save and out_file != '':
        plt.savefig(os.path.join(out_file, mode + '_result.png'), dpi = 500)
    plt.close()
    return
