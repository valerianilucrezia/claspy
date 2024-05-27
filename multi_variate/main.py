import sys
sys.path.append('../claspy')

import os
import argparse

from plot import *
from functions import *

from claspy.segmentation import BinaryClaSPSegmentation
from claspy.validation import map_validation_tests


distance = "euclidean_distance" 
window_size = 10 
n_jobs = 1

n_segments = BinaryClaSPSegmentation(distance=distance, window_size=window_size).n_segments
n_estimators = BinaryClaSPSegmentation(distance=distance, window_size=window_size).n_estimators
k_neighbours = BinaryClaSPSegmentation(distance=distance, window_size=window_size).k_neighbours
validation = BinaryClaSPSegmentation(distance=distance, window_size=window_size).validation
threshold = BinaryClaSPSegmentation(distance=distance, window_size=window_size).threshold
scored = BinaryClaSPSegmentation(distance=distance, window_size=window_size).score
early_stopping = BinaryClaSPSegmentation(distance=distance, window_size=window_size).early_stopping
excl_radius = BinaryClaSPSegmentation(distance=distance, window_size=window_size).excl_radius
random_state = BinaryClaSPSegmentation(distance=distance, window_size=window_size).random_state
validation_test = map_validation_tests(validation)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, default="/Users/lucreziavaleriani/Desktop/seq_res_snv_smooth.csv", help="input data")
    parser.add_argument('-out_dir', type=str, default="/Users/lucreziavaleriani/Desktop/", help="output dir")
    parser.add_argument('-mode', type=str, default="mult", help="type of cp choice bw: max, sum, mult")
    args = parser.parse_args()


    in_file = args.input
    output_dir = args.out_dir
    
    mode = args.mode
    vaf, baf, dr, bps = get_original_data_plot(in_file, 
                                               save = True,
                                               out_file = output_dir)
    bps = np.sort(bps)[1:-1]

    time_series = np.array([vaf, baf, dr])
    n_timepoints = time_series.shape[1]

    min_seg_size = window_size * excl_radius
    n_segments = time_series.shape[1] // min_seg_size

    dr_clasp, dr_profile, dr_cp, dr_range, dr_tree, dr_queue = get_first_cp(dr, n_estimators,
                    window_size, 
                    k_neighbours,
                    distance,
                    scored,
                    early_stopping,
                    excl_radius,
                    n_jobs, 
                    random_state, 
                    validation, 
                    threshold, 
                    n_segments)

    baf_clasp, baf_profile, baf_cp, baf_range, baf_tree, baf_queue  = get_first_cp(baf, n_estimators,
                    window_size, 
                    k_neighbours,
                    distance,
                    scored,
                    early_stopping,
                    excl_radius,
                    n_jobs, 
                    random_state, 
                    validation, 
                    threshold, 
                    n_segments)

    vaf_clasp, vaf_profile, vaf_cp, vaf_range, vaf_tree, vaf_queue = get_first_cp(vaf, n_estimators,
                    window_size, 
                    k_neighbours,
                    distance,
                    scored,
                    early_stopping,
                    excl_radius,
                    n_jobs, 
                    random_state, 
                    validation, 
                    threshold, 
                    n_segments)


    cp = take_first_cp(dr_profile, vaf_profile, baf_profile, mode)
    dr_tree, dr_queue, baf_tree, baf_queue, vaf_tree, vaf_queue, cp = validate_first_cp(cp, 
                                                                                        threshold, validation_test,
                                                                                        dr_clasp, dr_tree, dr_queue, dr_range, dr_profile,
                                                                                        baf_clasp, baf_tree, baf_queue, baf_range, baf_profile,
                                                                                        vaf_clasp, vaf_tree, vaf_queue, vaf_range, vaf_profile)


    CP =  find_cp_iterative(dr_clasp, dr_tree, dr_queue, dr_profile,
                        baf_clasp, baf_tree, baf_queue, baf_profile,
                        vaf_clasp, vaf_tree, vaf_queue, vaf_profile,
                        n_segments, validation, threshold, window_size, min_seg_size,
                        n_estimators, k_neighbours, distance, scored, early_stopping, 
                        excl_radius, n_jobs, random_state, n_timepoints, 
                        dr, baf, vaf, 
                        mode)

    plot_profile(dr, baf, vaf, bps, CP, 
                 title = mode, 
                 save = True, 
                 out_file = output_dir, 
                 mode = mode)