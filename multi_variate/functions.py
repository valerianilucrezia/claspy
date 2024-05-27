from queue import PriorityQueue
from claspy.utils import check_input_time_series, check_excl_radius
from claspy.clasp import ClaSPEnsemble
from claspy.validation import map_validation_tests
import numpy as np

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def get_first_cp(time_series, 
                 n_estimators,
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
                 n_segments):
    
    check_input_time_series(time_series)
    check_excl_radius(k_neighbours, excl_radius)

    n_timepoints = time_series.shape[0]
    min_seg_size = window_size * excl_radius
    
    
    queue = PriorityQueue()
    clasp_tree = []
    
    if n_segments == "learn":
        n_segments = time_series.shape[0] // min_seg_size
        
    prange = 0, time_series.shape[0]    
    clasp = ClaSPEnsemble(n_estimators=n_estimators,
                    window_size=window_size,
                    k_neighbours=k_neighbours,
                    distance=distance,
                    score=scored,
                    early_stopping=early_stopping,
                    excl_radius=excl_radius,
                    n_jobs=n_jobs,
                    random_state=random_state).fit(time_series, 
                                                   validation=validation, 
                                                   threshold=threshold)
            
    cp = my_split(clasp, clasp.profile, validation=validation, threshold=threshold) #clasp.split(validation=validation, threshold=threshold)
    # check if it is valid with significant test
    
    profile = clasp.profile
    return(clasp, profile, cp, prange, clasp_tree, queue)



def cp_is_valid(candidate, 
                change_points, 
                n_timepoints, 
                min_seg_size):
    value = [0] + change_points + [n_timepoints]
    
    for change_point in [0] + change_points + [n_timepoints]:
        left_begin = max(0, change_point - min_seg_size)
        right_end = min(n_timepoints, change_point + min_seg_size)
        
        if candidate in range(left_begin, right_end): 
            return False

    return True


def my_split(clasp, 
            profile, 
            sparse=True, 
            validation="significance_test", 
            threshold=1e-15):
    
    cp = np.argmax(profile)

    if validation is not None:
        validation_test = map_validation_tests(validation)
        if not validation_test(clasp, cp, threshold): 
            return 0 #None

    if sparse is True:
        return cp

    return time_series[:cp], time_series[cp:]

def local_segmentation(lbound, 
                       ubound, 
                       change_points, 
                       min_seg_size,
                       n_estimators,
                       window_size,
                       k_neighbours,
                       distance, 
                       score,
                       early_stopping,
                       excl_radius,
                       n_jobs,
                       random_states,
                       time_series,
                       validation,
                       threshold,
                       clasp_tree,
                       queue, 
                       n_timepoints):

    if ubound - lbound < 2 * min_seg_size: 
        return clasp_tree, queue, 0

    clasp = ClaSPEnsemble(
        n_estimators=n_estimators,
        window_size=window_size,
        k_neighbours=k_neighbours,
        distance=distance,
        score=score,
        early_stopping=early_stopping,
        excl_radius=excl_radius,
        n_jobs=n_jobs,
        random_state=random_states
    ).fit(time_series[lbound:ubound], validation=validation, threshold=threshold)

    cp = my_split(clasp, clasp.profile, validation=validation, threshold=threshold)
    if cp is None: 
        #print('cp is None')
        return clasp_tree, queue, 0
    
    score = clasp.profile[cp]

    if not cp_is_valid(lbound + cp, change_points, n_timepoints, min_seg_size):  #candidate, change_points, n_timepoints, min_seg_size
        #print('cp is not valid')
        return clasp_tree, queue, 0

    clasp_tree.append(((lbound, ubound), clasp))
    queue.put((-score, len(clasp_tree) - 1))    
    return clasp_tree, queue, score


def take_first_cp(dr_profile, 
                  vaf_profile, 
                  baf_profile, 
                  mode):    
    
    cps = [np.argmax(dr_profile), np.argmax(vaf_profile), np.argmax(baf_profile)]
    scores = [max(dr_profile), max(vaf_profile), max(baf_profile)]

    if mode == 'most_common':
        count = Counter(cps)  # Count the occurrences of each element
        most_common_cp, i = count.most_common(1)[0]  
        print('cp is:', most_common_cp)
        
    elif mode == 'max':
        scores = [max(dr_profile), max(vaf_profile), max(baf_profile)]
        idx = np.argmax(scores)
        most_common_cp = cps[idx]
        print('CP is:', most_common_cp)
        
        
    elif mode == 'mult':
        mult = dr_profile * vaf_profile * baf_profile
        most_common_cp = np.argmax(mult)
        print('CP is:', most_common_cp)
        
    elif mode == 'sum':
        sums = dr_profile + vaf_profile + baf_profile
        most_common_cp = np.argmax(sums/3)
        print('CP is:', most_common_cp)
        
    return most_common_cp
    
    
def validate_first_cp(cp,
                      threshold, validation_test,
                      dr_clasp, dr_tree, dr_queue, dr_range, dr_profile,
                      baf_clasp, baf_tree, baf_queue, baf_range, baf_profile,
                      vaf_clasp, vaf_tree, vaf_queue, vaf_range, vaf_profile):
    
    dr_val = validation_test(dr_clasp, cp, threshold)
    baf_val = validation_test(baf_clasp, cp, threshold)
    vaf_val = validation_test(vaf_clasp, cp, threshold)

    if dr_val or baf_val or vaf_val:
        print('CP passed:', cp)
        
        dr_tree.append((dr_range, dr_clasp))
        dr_queue.put((-dr_profile[cp], len(dr_tree) - 1))
        
        baf_tree.append((baf_range, baf_clasp))
        baf_queue.put((-baf_profile[cp], len(baf_tree) - 1))
        
        vaf_tree.append((vaf_range, vaf_clasp))
        vaf_queue.put((-vaf_profile[cp], len(vaf_tree) - 1))
        
        return dr_tree, dr_queue, baf_tree, baf_queue, vaf_tree, vaf_queue, cp
    
    else:
        return 
    
    
    
def find_cp_iterative(dr_clasp, dr_tree, dr_queue, dr_profile,
                      baf_clasp, baf_tree, baf_queue, baf_profile,
                      vaf_clasp, vaf_tree, vaf_queue, vaf_profile,
                      n_segments, validation, threshold, window_size, min_seg_size,
                      n_estimators, k_neighbours, distance, scored, early_stopping, 
                      excl_radius, n_jobs, random_state, n_timepoints, 
                      dr, baf, vaf, 
                      mode):
    CP = []
    dr_scores = []
    baf_scores = []
    vaf_scores = []
    
    for idx in range(n_segments - 1):
        if dr_queue.empty() and baf_queue.empty() and vaf_queue.empty(): 
            print('Stop because queue is empty')
            break
        
        if dr_queue.qsize() > 0:
            dr_priority, dr_clasp_tree_idx = dr_queue.get()
            (dr_lbound, dr_ubound), dr_clasp = dr_tree[dr_clasp_tree_idx]
            dr_cp = dr_lbound + my_split(dr_clasp, dr_clasp.profile, validation=validation, threshold=threshold)
            dr_profile[dr_lbound:dr_ubound - window_size + 1] = np.max([dr_profile[dr_lbound:dr_ubound - window_size + 1], dr_clasp.profile], axis=0)
            prof_dr = dr_clasp.profile
        
        else:
            dr_cp = 0
            dr_priority = 0
            prof_dr = np.array([])
            
        
        if baf_queue.qsize() > 0:
            baf_priority, baf_clasp_tree_idx = baf_queue.get()
            (baf_lbound, baf_ubound), baf_clasp = baf_tree[baf_clasp_tree_idx]
            baf_cp = baf_lbound + my_split(baf_clasp, baf_clasp.profile, validation=validation, threshold=threshold)
            baf_profile[baf_lbound:baf_ubound - window_size + 1] = np.max([baf_profile[baf_lbound:baf_ubound - window_size + 1], baf_clasp.profile], axis=0)
            prof_baf = baf_clasp.profile
        else:
            baf_cp = 0
            baf_priority = 0
            prof_baf = np.array([])
        
        
        if vaf_queue.qsize() > 0:
            vaf_priority, vaf_clasp_tree_idx = vaf_queue.get()
            (vaf_lbound, vaf_ubound), vaf_clasp = vaf_tree[vaf_clasp_tree_idx]
            vaf_cp = vaf_lbound + my_split(vaf_clasp, vaf_clasp.profile, validation=validation, threshold=threshold)
            vaf_profile[vaf_lbound:vaf_ubound - window_size + 1] = np.max([vaf_profile[vaf_lbound:vaf_ubound - window_size + 1], vaf_clasp.profile], axis=0)
            prof_vaf = vaf_clasp.profile
        else:
            vaf_cp = 0
            vaf_priority = 0
            prof_vaf = np.array([])
            
        
        
        all_bound = [(dr_lbound, dr_ubound), 
                     (baf_lbound, baf_ubound), 
                     (vaf_lbound, vaf_ubound)]
        all_cp = [dr_cp, baf_cp, vaf_cp]
        all_score = [dr_priority, baf_priority, vaf_priority]
        all_profile = [prof_dr, prof_baf, prof_vaf]
        all_profile = np.array([prof for prof in all_profile if prof.size != 0])
        
        
        if mode == 'most_common': 
            count = Counter(all_cp) 
            keep_cp, _ = count.most_common(1)[0]  
            
            
        elif mode == 'max': # max score
            idx = np.argmax(np.abs(all_score))
            keep_cp = all_cp[idx]
            
            
        elif mode == 'mult':
            new_score = np.prod(all_profile, axis = 0)#prof_dr * prof_baf * prof_vaf
            new_score[new_score == np.inf] = -10000
            keep_cp = dr_lbound + np.argmax(new_score)
            
            
        elif mode == 'sum':
            new_score = prof_dr + prof_baf + prof_vaf
            new_score[new_score == np.inf] = -10000
            keep_cp = dr_lbound + np.argmax(new_score/3)
            

        print(f'selected cp: {keep_cp}')
        CP.append(keep_cp)
        dr_scores.append(-dr_priority)
        baf_scores.append(-baf_priority)
        vaf_scores.append(-vaf_priority)
        
        # define new rrange
        lrange, rrange = (dr_lbound, keep_cp), (keep_cp, dr_ubound) 
        for prange in (lrange, rrange):
            
            low = prange[0]
            high = prange[1]
            dr_tree, dr_queue, dr_scores_tmp = local_segmentation(low, 
                                                                high, 
                                                                CP, 
                                                                min_seg_size,
                                                                n_estimators,
                                                                window_size,
                                                                k_neighbours,
                                                                distance, 
                                                                scored,
                                                                early_stopping,
                                                                excl_radius,
                                                                n_jobs,
                                                                random_state,
                                                                dr, 
                                                                validation, 
                                                                threshold,
                                                                dr_tree, 
                                                                dr_queue,
                                                                n_timepoints)
            
            baf_tree, baf_queue, baf_scores_tmp = local_segmentation(low, 
                                                                    high, 
                                                                    CP, 
                                                                    min_seg_size,
                                                                    n_estimators,
                                                                    window_size,
                                                                    k_neighbours,
                                                                    distance, 
                                                                    scored,
                                                                    early_stopping,
                                                                    excl_radius,
                                                                    n_jobs,
                                                                    random_state,
                                                                    baf, 
                                                                    validation, 
                                                                    threshold,
                                                                    baf_tree, 
                                                                    baf_queue,
                                                                    n_timepoints)
            
            vaf_tree, vaf_queue, vaf_scores_tmp = local_segmentation(low, 
                                                                    high, 
                                                                    CP, 
                                                                    min_seg_size,
                                                                    n_estimators,
                                                                    window_size,
                                                                    k_neighbours,
                                                                    distance, 
                                                                    scored,
                                                                    early_stopping,
                                                                    excl_radius,
                                                                    n_jobs,
                                                                    random_state,
                                                                    vaf, 
                                                                    validation, 
                                                                    threshold,
                                                                    vaf_tree, 
                                                                    vaf_queue,
                                                                    n_timepoints)
            
            #print(f'scores = {dr_scores_tmp}, {baf_scores_tmp}, {vaf_scores_tmp}')
            if dr_scores_tmp != None:
                dr_scores.append(dr_scores_tmp)
                
            if baf_scores_tmp != None:
                baf_scores.append(baf_scores_tmp)
                
            if vaf_scores_tmp != None:
                vaf_scores.append(vaf_scores_tmp)
    
    return CP
    