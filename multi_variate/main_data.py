import sys
sys.path.append('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy')

import os
import numpy as np
import argparse

from multi_variate.data_import import *
from multi_variate.multivariate_segmentation import *

from claspy.segmentation import BinaryClaSPSegmentation
from claspy.validation import map_validation_tests

thr = [1e-15, 1e-10, 1e-5]
wsize = [5, 10]


distance = "euclidean_distance" 
n_jobs = 1

n_segments = BinaryClaSPSegmentation(distance=distance).n_segments
n_estimators = BinaryClaSPSegmentation(distance=distance).n_estimators
k_neighbours = BinaryClaSPSegmentation(distance=distance).k_neighbours
validation = BinaryClaSPSegmentation(distance=distance).validation
scored = BinaryClaSPSegmentation(distance=distance).score
early_stopping = BinaryClaSPSegmentation(distance=distance).early_stopping
excl_radius = BinaryClaSPSegmentation(distance=distance).excl_radius
random_state = BinaryClaSPSegmentation(distance=distance).random_state
validation_test = map_validation_tests(validation)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type = str, default ="/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_nanopore/data/smooth_data.tsv", help = "input path")
    parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
    parser.add_argument('-output', type = str, default = "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_nanopore/claspy/", help = "output path")
    
    args = parser.parse_args() 
    mode = args.mode
    in_file = args.input
    out_dir = args.output
        
    print('Mode:', mode)
    output_dir = out_dir
    os.makedirs(output_dir, exist_ok=True)
    
    for window_size in wsize:
        print('Window size:', window_size)
        
        for threshold in thr:
            print('Threshold:', threshold)
            name = f'{mode}_{window_size}_{threshold}'
            vaf, baf, dr, maf = get_data(in_file)
            
            
            time_series = np.array([vaf, baf, dr])
            n_timepoints = time_series.shape[1]

            min_seg_size = window_size * excl_radius
            n_segments = time_series.shape[1] // min_seg_size

            try:
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
                
                if type(maf) == type(np.array([1,2,3])):
                    maf_clasp, maf_profile, maf_cp, maf_range, maf_tree, maf_queue = get_first_cp(maf, n_estimators,
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
                    
                    cp = take_first_cp(dr_profile, vaf_profile, baf_profile, mode, maf_profile)
                    dr_tree, dr_queue, baf_tree, baf_queue, vaf_tree, vaf_queue, maf_tree, maf_queue, cp = validate_first_cp(cp, threshold, validation_test,
                                                                                                                dr_clasp, dr_tree, dr_queue, dr_range, dr_profile,
                                                                                                                baf_clasp, baf_tree, baf_queue, baf_range, baf_profile,
                                                                                                                vaf_clasp, vaf_tree, vaf_queue, vaf_range, vaf_profile,
                                                                                                                maf_clasp, maf_tree, maf_queue, maf_range, maf_profile)

                    
                    CP =  MAF_find_cp_iterative(dr_clasp, dr_tree, dr_queue, dr_profile,
                                            baf_clasp, baf_tree, baf_queue, baf_profile,
                                            vaf_clasp, vaf_tree, vaf_queue, vaf_profile,
                                            n_segments, validation, threshold, window_size, min_seg_size,
                                            n_estimators, k_neighbours, distance, scored, early_stopping, 
                                            excl_radius, n_jobs, random_state, n_timepoints, 
                                            dr, baf, vaf,
                                            mode,
                                            maf_clasp, maf_tree, maf_queue, maf_profile, maf)
                    
                    np.save(file = f'{output_dir}/{name}_maf.npy', arr = np.array(CP))
                else:

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
                        
                    np.save(file = f'{output_dir}/{name}.npy', arr = np.array(CP))
            except:
               print(f'Not pass: {name}')
               pass
