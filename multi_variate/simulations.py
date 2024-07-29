from multivariate_clasp import MultivariateClaSP
import argparse
import sys
import os
import numpy as np

sys_path = None

def check_req_args(input, mode, sys_path, sim):
    if os.path.isdir(input) is False:
        raise TypeError
    if mode not in ["max", "sum", "mult"]:
        raise ValueError
    if sim is None:
        raise TypeError

def run_simulations(input, mode, output, sys_path, sim, frequencies, n_segments, n_estimators, window_size, k_neighbors,
                      distance, score, early_stopping, validation, threshold, excl_radius, n_jobs, random_state):
    if sys_path is not None:
        sys.path.append(sys_path)
        

    
    # FIXME sim not necessary here. can still be used to name output
    base = os.path.join(input, sim)
    combinations = os.listdir(base)

    for c in combinations:
        path = os.path.join(base, c)
        if os.path.isdir(path) and c[0] == 'c':
            print('Running:', c)

            in_file = os.path.join(path, 'smooth_snv.csv')   
            output_dir = os.path.join(output, sim, c)
            os.makedirs(output_dir, exist_ok=True)
            
            for wsize in window_size:
                
                for thr in threshold:

                    multivariate_ts = MultivariateClaSP(in_file, mode, output_dir, frequencies, n_segments, n_estimators, wsize, k_neighbors,
                      distance, score, early_stopping, validation, thr, excl_radius, n_jobs, random_state)
                    multivariate_ts.analyze_time_series()

                    # save out results
                    np.save(file = os.path.join(multivariate_ts.out_dir, f'{multivariate_ts.name}_bps.npy'), arr = np.array(multivariate_ts.bps))
                    np.save(file = os.path.join(multivariate_ts.out_dir, f'{multivariate_ts.name}_cp.npy'), arr = np.array(multivariate_ts.CP))

                    multivariate_ts.plot_profile(save=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type = str, help = "input path to simulation data folders")
    parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
    # optional args to be passed
    parser.add_argument('--output', type = str, default=None, help = "output path")
    parser.add_argument('--sys_path', type=str, default=None, help="Path to claspy on HPC system")
    parser.add_argument('--sim', type=str, default="Sim", help="Simulation run name")
    parser.add_argument('--frequencies', type=list, default=["vaf", "baf", "dr"])
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

    run_simulations(args.input, args.mode, args.output, args.sys_path, args.sim, args.frequencies, args.n_segments, args.n_estimators, args.window_size, args.k_neighbors,
                    args.distance, args.score, args.early_stopping, args.validation, args.threshold, args.excl_radius, args.n_jobs, args.random_state)