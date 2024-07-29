import sys
# FIXME if sys.path.append is important should move to be function argument
# sys.path.append('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy')

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from data_import import get_data_csv, get_data_tsv
from multivariate_segmentation import MultivariateClaSPSegmentation, take_first_cp, validate_first_cp, find_cp_iterative

# TODO inherit from some claspy class to run error checking on kwargs? Is argparse type enough
class MultivariateClaSP:
    def __init__(self, input, mode, out_dir, frequencies=["vaf", "baf", "dr"], n_segments="learn", n_estimators=10, window_size=5, k_neighbours=3, distance="euclidean_distance", score="roc_auc", early_stopping=True, validation="significance_test", threshold=1e-15, excl_radius=5, n_jobs=1, random_state=2357):
        # take dict of inherited args to pass to MultivariateClaSPSegmentation. Delete params only used with this class
        self.kwargs = locals()
        del self.kwargs["self", "input", "mode", "out_dir", "frequencies"]
        # assign as usual
        self.frequencies = frequencies
        self.input = input
        self.mode = mode
        self.out_dir = out_dir

        self.get_data = None

        self._check_args()

    def _check_args(self):
        # TODO more specific error checking
        if os.path.isfile(self.input) is True:
            if self.input.endswith(".csv"):
                self.get_data = get_data_csv
            elif self.input.endswith(".tsv"):
                self.get_data = get_data_tsv
        else:
            raise TypeError(f"input file must be csv or tsv, not {self.input}")
        # should be all that is needed as input is already checked
        if self.out_dir is None:
            self.out_dir = os.path.dirname(self.input)
        if type(self.mode) == str or self.mode in ["max", "sum", "mult"] is False:
            raise TypeError(f"mode must be string of one of the following options: max, sum, mult")
        # if all(self.kwargs["frequencies"], str) is False:
        #     raise TypeError(f"If frequencies is specified, list items must be strings")


    def analyze_time_series(self):
        self.name = f'{self.mode}_{self.window_size}_{self.threshold}'

        # returns dictionary to preserve variable names
        original_data = self.get_data(self.input)

        # sort bps
        self.bps = np.sort(original_data["bps"])[1:-1]

        time_series = np.array([i for i in original_data[i] if i in self.frequencies])
        # n_timepoints = time_series.shape[1]

        # now done in multivariateClaSP init
        # min_seg_size = window_size * excl_radius
        # n_segments = time_series.shape[1] // min_seg_size
        try:
            self.multivariate_clasp_objects = {}
            for i in self.frequencies:
                ts_obj = MultivariateClaSPSegmentation(time_series=original_data[i], **self.kwargs)
                ts_obj.n_timepoints = time_series.shape[1]
                ts_obj.n_segments = ts_obj.n_timepoints // ts_obj.min_seg_size
                ts_obj.get_first_cp()
                self.multivariate_clasp_objects[i] = ts_obj


            # call with profiles for each time series
            cp = take_first_cp(self.multivariate_clasp_objects,
                                self.mode)
            # feed in created objects to be checked
            validate_first_cp(multivariate_clasp_objects=self.multivariate_clasp_objects, cp=cp)
            
            self.CP = find_cp_iterative(self.multivariate_clasp_objects, self.mode)

        except:
            print(f'Not passed: {self.name}')
            pass
    
    @staticmethod
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

    def get_original_data_plot(self, save: bool = False):
        snv = pd.read_csv(self.input)        
        self.plot_data(snv)       
        if save:
            plt.savefig(fname = os.path.join(self.out_dir, 'original_data.png'))
        plt.close()


    def plot_profile(self, pred_bps = [],
                    title = '', 
                    save = False):
        
        # fig = plt.figure(figsize=(10, 6))
        fig, axs = plt.subplots(nrows = 3, ncols = 1, figsize=(10, 8))
        # ax = plt.subplot(111)
        variables = self.multivariate_clasp_objects.keys()
        for i in range(0, len(variables)-1):
            profile = self.multivariate_clasp_objects[variables[i]].profile
            axs[i].plot(profile, label = variables)
            axs[i].set_ylabel(variables[i])
            axs[i].vlines(self.bps, ymin = min(profile) - 0.05, ymax = max(profile) + 0.05, colors = 'tab:green', label = 'True BP', linestyles = 'dashed')
            axs[i].vlines(pred_bps, ymin = min(profile), ymax = max(profile), colors = 'tab:olive',  label = 'Predicted BP')


        #axs.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
        fig.suptitle(title)
        fig.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.out_dir, f'{self.mode}_{title}_result.png'), dpi = 500)
        plt.close()


# TODO probably best to move this to separate entrance module now
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-input', type = str, help = "input path to smoothed data csv file")
#     parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
#     parser.add_argument('-output', type = str, default=None, help = "output path")
#     # optional args to be passed to ClaSP, called in main as **kwargs
#     parser.add_argument('--frequencies', type=list, default=["vaf", "baf", "dr"])
#     parser.add_argument('--n_segments', type = str, default = "learn")
#     parser.add_argument('--n_estimators', type = int, default = 10)
#     parser.add_argument('--window_size', type = int, default = 5)
#     parser.add_argument('--k_neighbors', type = int, default = 3)
#     parser.add_argument('--distance', type = str, default =  "euclidean_distance")
#     parser.add_argument('--score', type = str, default = "roc_auc")
#     parser.add_argument('--early_stopping', type = bool, default = True)
#     parser.add_argument('--validation', type = str, default = "significance_test")
#     parser.add_argument('--threshold', type = float, default = 1e-15)
#     parser.add_argument('--excl_radius', type = int, default = 5)
#     parser.add_argument('--n_jobs', type = int, default = 1)
#     parser.add_argument('--random_state', type = int, default = 2357)
    
#     args = parser.parse_args()
#     # gross
#     MultivariateClaSP(args.input, args.mode, args.output, args.frequencies, args.n_segments, args.n_estimators, args.window_size, args.k_neighbors,
#                       args.distance, args.score, args.early_stopping, args.validation, args.threshold, args.excl_radius, args.n_jobs, args.random_state).analyze_time_series()