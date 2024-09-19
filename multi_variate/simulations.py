from multivariate_clasp import MultivariateClaSP
import argparse
import sys
import os
import numpy as np

sys_path = None

def check_req_args(input, mode, sys_path, sim) -> None:
    """_summary_

    Parameters
    ----------
    input : str
        _description_
    mode : _type_
        _description_
    sys_path : _type_
        _description_
    sim : _type_
        _description_

    Raises
    ------
    TypeError
        _description_
    ValueError
        _description_
    TypeError
        _description_
    """
    if os.path.isdir(input) is False:
        raise TypeError
    if mode not in ["max", "sum", "mult"]:
        raise ValueError
    if sim is None:
        raise TypeError

def run_simulations(input: str, mode: str, output: str, sys_path: str, sim: str, frequencies: list[str], n_segments: str | int, n_estimators: int, window_size: str | int,
                    k_neighbors: int, distance: str, score: str, early_stopping: bool, validation: str, threshold: float, excl_radius: int, n_jobs: int, random_state: int):
    """Run simulations given a file of allele frequencies and multiple window_size and threshold values to check and save out a plot with identified change points.

    Parameters
    ----------
    input : str
        The path to a csv file holding gene frequencies (vaf, baf, dr, etc.).

    mode : str
        Denotes the method for deriving a change point, either sum, max, or multi.
    
    output : str
        The path to store outputs.

    sys_path : str
        The path to store outputs.
    
    sim : str
        Simulation run name.

    frequencies : list
        Strings denoting the gene frequencies to be used. default is ["vaf", "median_baf", "median_dr"]

    n_segments : str or int, default="learn"
        The number of segments to split the time series into. By default, the numbers
        of segments is inferred automatically by applying a change point validation test.

    n_estimators : int, default=10
        The number of ClaSPs in the ensemble.

    window_size : str or int, default="suss"
        The window size detection method or size of the sliding window used in the ClaSP
        algorithm. Valid implementations include: 'suss', 'fft', and 'acf'.

    k_neighbours : int, default=3
        The number of nearest neighbors to use in the ClaSP algorithm.

    distance: str, default="euclidean_distance"
        The name of the distance function to be computed for determining the k-NNs. Available
        options are "znormed_euclidean_distance" and "euclidean_distance".

    score : str, default="roc_auc"
        The name of the scoring metric to use in ClaSP. Available options are "roc_auc",
        "f1".

    early_stopping : bool, optional
        Determines if ensembling is stopped, once a validated change point is found or
        the ClaSP models do not improve anymore. Default is True.

    validation : str, optional
        The validation method to use for determining the significance of the change point.
        The available methods are "significance_test" and "score_threshold". Default is
        "significance_test".

    threshold : float, optional
        The threshold value to use for the validation test. If the validation method is
        "significance_test", this value represents the p-value threshold for rejecting the
        null hypothesis. If the validation method is "score_threshold", this value represents
        the threshold score for accepting the change point. Default is 1e-15.

    excl_radius : int, default=5
        The radius (in multiples of the window size) around each point in the time series to exclude
        when searching for change points.

    n_jobs : int, optional (default=1)
        Amount of threads used in the ClaSP computation.
    
    random_state : int, optional
        Sets random seed for reproducibility. Default is 2357.
    """
    
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

                    try:

                        multivariate_ts = MultivariateClaSP(in_file, mode, output_dir, frequencies, n_segments, n_estimators, wsize, k_neighbors,
                        distance, score, early_stopping, validation, thr, excl_radius, n_jobs, random_state)
                        multivariate_ts.analyze_time_series()

                        # save out results
                        np.save(file = os.path.join(multivariate_ts.out_dir, f'{multivariate_ts.name}_bps.npy'), arr = np.array(multivariate_ts.bps))
                        np.save(file = os.path.join(multivariate_ts.out_dir, f'{multivariate_ts.name}_cp.npy'), arr = np.array(multivariate_ts.CP))

                        multivariate_ts.plot_profile(save=True)
                    except:
                        print(f"{wsize}_{thr} not passed")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type = str, help = "input path to simulation data folders")
    parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
    # optional args to be passed
    parser.add_argument('--output', type = str, default=None, help = "output path")
    parser.add_argument('--sys_path', type=str, default=None, help="Path to claspy on HPC system")
    parser.add_argument('--sim', type=str, default="Sim", help="Simulation run name")
    parser.add_argument('--frequencies', type=list, default=["vaf", "median_baf", "median_dr"])
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