import sys
sys.path.append('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/claspy')

import os
import numpy as np
import argparse

from plot import *
from functions import *

from claspy.segmentation import BinaryClaSPSegmentation
from claspy.validation import map_validation_tests

thr = [1e-15, 1e-10, 1e-5]
wsize = [5,10,50]
#threshold = 1e-15
#window_size = 10 


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
    parser.add_argument('-input', type = str, default ="/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/sim_data_races/data_races/", help = "input path")
    parser.add_argument('-mode', type = str, default = "mult", help= "type of cp choice bw: max, sum, mult")
    parser.add_argument('-sim', type = str, default = "sim_3", help= "simulation name")
    parser.add_argument('-output', type = str, default = "/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_races/claspy/", help = "output path")
    
    args = parser.parse_args() 
    mode = args.mode
    sim = args.sim
    in_dir = args.input
    out_dir = args.output
    
    base = os.path.join(in_dir, sim)
    combinations = os.listdir(base) 
    
    print('Simulation:', sim)
    print('Mode:', mode)
    
    for c in combinations:
        path = os.path.join(base, c)
        if os.path.isdir(path) and c[0] == 'c':
            print('Running:', c)

            in_file = os.path.join(path, 'smooth_snv.csv')   
            output_dir = os.path.join(out_dir, sim, c)
            os.makedirs(output_dir, exist_ok=True)
            
            for window_size in wsize:
                print('Window size:', window_size)
                
                for threshold in thr:
                    print('Threshold:', threshold)
    
                    name = f'{mode}_{window_size}_{threshold}'
                    vaf, baf, dr, bps = get_original_data_plot(in_file, 
                                                            save = True,
                                                            out_file = output_dir)
                    bps = np.sort(bps)[1:-1]
                    np.save(file = f'{os.path.join(out_dir, sim)}/bp.npy', arr = np.array(bps))

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
                        plot_profile(dr, baf, vaf, bps, CP, 
                                            title = f'{c}-{name}', 
                                            save = True, 
                                            out_file = output_dir, 
                                            mode = f'{mode}_{name}')
                    except:
                        print(f'Not passe: {name}')
                        pass
