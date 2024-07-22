import sys
sys.path.append('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy')

import os
import numpy as np
import argparse

from plot import *
from functions import *

# TODO: Why are these constants? Why is the main function looping through
# thr = [1e-15, 1e-10, 1e-5]
# wsize = [5, 10, 50]
#threshold = 1e-15
#window_size = 10

def main(in_dir, mode, sim, out_dir, **kwargs):
    base = os.path.join(in_dir, sim)
    # TODO: what is combinations?
    combinations = os.listdir(base) 
    
    print('Simulation:', sim)
    print('Mode:', mode)

    # take threshold and window size args from kwargs and keep as list, then delete from kwargs to avoid overwriting
    thr = kwargs['threshold']
    wsize = kwargs['window_size']
    del kwargs['threshold']
    del kwargs['window_size'] 
    
    for c in combinations:
        path = os.path.join(base, c)
        if os.path.isdir(path) and c[0] == 'c':
            print('Running:', c)
            # FIXME: add functionality for custom filename? or enforce smooth_snv.csv naming?
            in_file = os.path.join(path, 'smooth_snv.csv')   
            output_dir = os.path.join(out_dir, sim, c)
            os.makedirs(output_dir, exist_ok=True)
            # FIXME: kwargs will overwite window_size... but I can't ovverwrite kwargs["window_size"]...
            for window_size in wsize:
                print('Window size:', window_size)
                
                for threshold in thr:
                    print('Threshold:', threshold)
    
                    name = f'{mode}_{window_size}_{threshold}'
                    # returns dictionary now to preserve variable names
                    original_data = get_original_data_plot(in_file, 
                                                            save = True,
                                                            out_file = output_dir)
                    bps = np.sort(original_data[-1])[1:-1]
                    print(bps)
                    np.save(file = f'{os.path.join(out_dir, sim)}/{c}_bp.npy', arr = np.array(bps))
                    # TODO: change to acess via original data list
                    time_series = np.array([i for i in original_data.values()[:-1]])
                    # n_timepoints = time_series.shape[1]

                    # now done in multivariateClaSP init
                    # min_seg_size = window_size * excl_radius
                    # n_segments = time_series.shape[1] // min_seg_size
                    
                    try:
                        multivariate_clasp_objects = {}
                        for i in range(0, len(original_data)-2):
                            ts_obj = multivariateClaSP(original_data[original_data.keys()[i]], window_size=window_size, threshold=threshold, **kwargs)
                            ts_obj.n_timepoints = time_series.shape[1]
                            ts_obj.n_segments = ts_obj.n_timepoints // ts_obj.min_seg_size
                            ts_obj.get_first_cp()
                            multivariate_clasp_objects[original_data.keys()[i]] = ts_obj


                        
                        # call with profiles for each time series
                        # TODO: is there a better way to do this?? parse dict values into variables of name key
                        cp = take_first_cp(multivariate_clasp_objects["dr"].profile,
                                           multivariate_clasp_objects["baf"].profile,
                                           multivariate_clasp_objects["vaf"].profile,
                                           mode)
                        # feed in created objects to be checked
                        validate_first_cp(multivariate_clasp_objects=multivariate_clasp_objects.values(), cp=cp)
                        
                        CP = find_cp_iterative(multivariate_clasp_objects, mode)
                            
                        np.save(file = f'{output_dir}/{name}.npy', arr = np.array(CP))
                        plot_profile(multivariate_clasp_objects, bps, CP, 
                                            title = f'{c}-{name}', 
                                            save = True, 
                                            out_file = output_dir, 
                                            mode = f'{mode}_{name}')
                    except:
                        print(f'Not passe: {name}')
                        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type = str, default ="/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/", help = "input path")
    parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
    parser.add_argument('-sim', type = str, default = "sim_3", help= "simulation name")
    parser.add_argument('-output', type = str, default = "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_races/claspy/", help = "output path")
    # optional args to be passed to ClaSP, called in main as **kwargs
    # TODO: confirm list type for window size and threshold
    parser.add_argument('--n_segments', type = str, default = "learn")
    parser.add_argument('--n_estimators', type = int, default = 10)
    parser.add_argument('--window_size', type = list, default = [5, 10, 50])
    parser.add_argument('--k_neighbors', type = int, default = 3)
    parser.add_argument('--distance', type = str, default =  "euclidean_distance")
    parser.add_argument('--score', type = str, default = "roc_auc")
    parser.add_argument('--early_stopping', type = bool, default = True)
    parser.add_argument('--validation', type = str, default = "significance_test")
    parser.add_argument('--threshold', type = list, default = [1e-15, 1e-10, 1e-5])
    parser.add_argument('--excl_radius', type = int, default = 5)
    parser.add_argument('--n_jobs', type = int, default = 1)
    parser.add_argument('--random_state', type = int, default = 2357)
    
    args = parser.parse_args() 
    main(args.input, args.mode, args.sim, args.output, args.n_segments, args.n_estimators, args.window_size, args.k_neighbors, args.distance, args.score, args.early_stopping,
         args.validation, args.threshold, args.excl_radius, args.n_jobs, args.random_state)