from queue import PriorityQueue
from claspy.utils import check_input_time_series, check_excl_radius
from claspy.clasp import ClaSPEnsemble
from claspy.validation import map_validation_tests
from claspy.validation import significance_test_2
from claspy.validation import significance_test
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
    
    
    queue = [] #PriorityQueue()
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
    # print('first cp',cp)
    # # check if it is valid with significant test
    
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
   # print(cp)

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
        return clasp_tree, queue, 0
    
    score = clasp.profile[cp]
    #print('cp, score=', cp, score)
    
    if not cp_is_valid(lbound + cp, change_points, n_timepoints, min_seg_size):  #candidate, change_points, n_timepoints, min_seg_size
        #print('cp is not valid')
        #return clasp_tree, queue, 0
        score = 0

    clasp_tree.append(((lbound, ubound), clasp))
    #queue.put((-score, len(clasp_tree) - 1))    
    queue.append((-score, len(clasp_tree) - 1))
    return clasp_tree, queue, score


def take_first_cp(dr_profile, 
                  vaf_profile, 
                  baf_profile, 
                  mode, 
                  maf_profile = []):    
    
    if maf_profile != []:
        cps = [np.argmax(dr_profile), np.argmax(vaf_profile), np.argmax(baf_profile), np.argmax(maf_profile)]
        scores = [max(dr_profile), max(vaf_profile), max(baf_profile), max(maf_profile)]
        
        if mode == 'max':
            idx = np.argmax(scores)
            most_common_cp = cps[idx]
            
        elif mode == 'mult':
            mult = dr_profile * vaf_profile * baf_profile * maf_profile
            most_common_cp = np.argmax(mult)
            
        elif mode == 'sum':
            sums = dr_profile + vaf_profile + baf_profile + maf_profile
            most_common_cp = np.argmax(sums/4)
        
        
    else:
        cps = [np.argmax(dr_profile), np.argmax(vaf_profile), np.argmax(baf_profile)]
        scores = [max(dr_profile), max(vaf_profile), max(baf_profile)]
            
        if mode == 'max':
            idx = np.argmax(scores)
            most_common_cp = cps[idx]
            
        elif mode == 'mult':
            mult = dr_profile * vaf_profile * baf_profile
            most_common_cp = np.argmax(mult)
            
        elif mode == 'sum':
            sums = dr_profile + vaf_profile + baf_profile
            most_common_cp = np.argmax(sums/3)
        
    return most_common_cp
    
    
def validate_first_cp(cp,
                      threshold, validation_test,
                      dr_clasp, dr_tree, dr_queue, dr_range, dr_profile,
                      baf_clasp, baf_tree, baf_queue, baf_range, baf_profile,
                      vaf_clasp, vaf_tree, vaf_queue, vaf_range, vaf_profile,
                      maf_clasp = [], maf_tree = [], 
                      maf_queue = [], maf_range = [], 
                      maf_profile = []):
    
    dr_val = validation_test(dr_clasp, cp, threshold)
    baf_val = validation_test(baf_clasp, cp, threshold)
    vaf_val = validation_test(vaf_clasp, cp, threshold)
    
    if maf_clasp != []:
        maf_val = validation_test(maf_clasp, cp, threshold)
        
        if dr_val or baf_val or vaf_val or maf_val: 
            dr_tree.append((dr_range, dr_clasp))
            dr_queue.append((-dr_profile[cp], len(dr_tree) - 1))
            
            baf_tree.append((baf_range, baf_clasp))
            baf_queue.append((-baf_profile[cp], len(baf_tree) - 1))
            
            vaf_tree.append((vaf_range, vaf_clasp))
            vaf_queue.append((-vaf_profile[cp], len(vaf_tree) - 1))
            
            maf_tree.append((maf_range, maf_clasp))
            maf_queue.append((-maf_profile[cp], len(maf_tree) - 1))
            
            return dr_tree, dr_queue, baf_tree, baf_queue, vaf_tree, vaf_queue, maf_tree, maf_queue, cp
        
        else:
            return
        
    else:
    
        if dr_val or baf_val or vaf_val:        
            dr_tree.append((dr_range, dr_clasp))
            dr_queue.append((-dr_profile[cp], len(dr_tree) - 1))
            
            baf_tree.append((baf_range, baf_clasp))
            baf_queue.append((-baf_profile[cp], len(baf_tree) - 1))
            
            vaf_tree.append((vaf_range, vaf_clasp))
            vaf_queue.append((-vaf_profile[cp], len(vaf_tree) - 1))
            
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
                      mode,
                      maf_clasp = [], maf_tree = [], maf_queue = [], maf_profile = [], maf = []):
    CP = []
    dr_scores = []
    baf_scores = []
    vaf_scores = []
    
    for idx in range(n_segments - 1):
        if (len(dr_queue) == 0 and len(baf_queue) == 0) and len(vaf_queue) == 0:
            print('Stop because queue is empty')
            break
        
        if len(dr_queue) > 0:
            dr_priority, dr_clasp_tree_idx = dr_queue.pop() #dr_queue.get()
            (dr_lbound, dr_ubound), dr_clasp = dr_tree[dr_clasp_tree_idx]
            dr_cp = dr_lbound + my_split(dr_clasp, dr_clasp.profile, validation=validation, threshold=threshold)
            dr_profile[dr_lbound:dr_ubound - window_size + 1] = np.max([dr_profile[dr_lbound:dr_ubound - window_size + 1], dr_clasp.profile], axis=0)
            prof_dr = dr_clasp.profile
        
        else:
            dr_cp = 0
            dr_priority = 0
            prof_dr = np.array([])
        

        #if baf_queue.qsize() > 0:
        if len(baf_queue) > 0: 
            baf_priority, baf_clasp_tree_idx = baf_queue.pop() #baf_queue.get()
            (baf_lbound, baf_ubound), baf_clasp = baf_tree[baf_clasp_tree_idx]
            baf_cp = baf_lbound + my_split(baf_clasp, baf_clasp.profile, validation=validation, threshold=threshold)
            baf_profile[baf_lbound:baf_ubound - window_size + 1] = np.max([baf_profile[baf_lbound:baf_ubound - window_size + 1], baf_clasp.profile], axis=0)
            prof_baf = baf_clasp.profile
        
        else:
            baf_cp = 0
            baf_priority = 0
            prof_baf = np.array([])
        
        #if vaf_queue.qsize() > 0:
        if len(vaf_queue) > 0:
            vaf_priority, vaf_clasp_tree_idx = vaf_queue.pop() #vaf_queue.get()
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
        
                    
        if mode == 'max': # max score
            idx = np.argmax(np.abs(all_score))
            keep_cp = all_cp[idx]
            
            
        elif mode == 'mult':
            new_score = np.prod(all_profile, axis = 0) #prof_dr * prof_baf * prof_vaf
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score)
            keep_cp = dr_lbound + np.argmax(new_score)
            
            
        elif mode == 'sum':
            new_score = np.sum(all_profile, axis = 0)
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score/3)
            keep_cp = dr_lbound + np.argmax(new_score/3)
            
        
        if mode != 'max':
            all_val = [] 
            if prof_vaf.size > 0:
                val_vaf = significance_test(vaf_clasp, test_cp, threshold)
                all_val.append(all_val)
            
            if prof_baf.size > 0:
                val_baf = significance_test(baf_clasp, test_cp, threshold)
                all_val.append(val_baf)
            
            if prof_dr.size > 0:
                val_dr = significance_test(dr_clasp, test_cp, threshold)
                all_val.append(val_dr)
            
            if (True in all_val):
                CP.append(keep_cp)
        
        else:
            CP.append(keep_cp)
            
        
        
        dr_scores.append(-dr_priority)
        baf_scores.append(-baf_priority)
        vaf_scores.append(-vaf_priority)
        
        # define new rrange
        lrange, rrange = (dr_lbound, keep_cp), (keep_cp, dr_ubound) 
        for prange in (lrange, rrange):
            #print('\n Search in bound:', prange)
            
            low = prange[0]
            high = prange[1]

            #print('\n DR')
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
            #print('\n BAF')
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

            #print('\n VAF')
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
            
            # print(f'scores = {dr_scores_tmp}, {baf_scores_tmp}, {vaf_scores_tmp}')
            if dr_scores_tmp != None:
                dr_scores.append(dr_scores_tmp)
                
            if baf_scores_tmp != None:
                baf_scores.append(baf_scores_tmp)
                
            if vaf_scores_tmp != None:
                vaf_scores.append(vaf_scores_tmp)
    
    return CP

def MAF_find_cp_iterative(dr_clasp, dr_tree, dr_queue, dr_profile,
                      baf_clasp, baf_tree, baf_queue, baf_profile,
                      vaf_clasp, vaf_tree, vaf_queue, vaf_profile,
                      n_segments, validation, threshold, window_size, min_seg_size,
                      n_estimators, k_neighbours, distance, scored, early_stopping, 
                      excl_radius, n_jobs, random_state, n_timepoints, 
                      dr, baf, vaf, 
                      mode,
                      maf_clasp, maf_tree, maf_queue, maf_profile, maf):
    CP = []
    dr_scores = []
    baf_scores = []
    vaf_scores = []
    maf_scores = []
    
    for idx in range(n_segments - 1):
        if len(dr_queue) == 0 and len(baf_queue) == 0 and len(vaf_queue) == 0 and len(maf_queue) == 0:
            print('Stop because queue is empty')
            break
        
        if len(dr_queue) > 0:
            dr_priority, dr_clasp_tree_idx = dr_queue.pop() #dr_queue.get()
            (dr_lbound, dr_ubound), dr_clasp = dr_tree[dr_clasp_tree_idx]
            dr_cp = dr_lbound + my_split(dr_clasp, dr_clasp.profile, validation=validation, threshold=threshold)
            dr_profile[dr_lbound:dr_ubound - window_size + 1] = np.max([dr_profile[dr_lbound:dr_ubound - window_size + 1], dr_clasp.profile], axis=0)
            prof_dr = dr_clasp.profile
        
        else:
            dr_cp = 0
            dr_priority = 0
            prof_dr = np.array([])
        

        #if baf_queue.qsize() > 0:
        if len(baf_queue) > 0: 
            baf_priority, baf_clasp_tree_idx = baf_queue.pop() #baf_queue.get()
            (baf_lbound, baf_ubound), baf_clasp = baf_tree[baf_clasp_tree_idx]
            baf_cp = baf_lbound + my_split(baf_clasp, baf_clasp.profile, validation=validation, threshold=threshold)
            baf_profile[baf_lbound:baf_ubound - window_size + 1] = np.max([baf_profile[baf_lbound:baf_ubound - window_size + 1], baf_clasp.profile], axis=0)
            prof_baf = baf_clasp.profile
        
        else:
            baf_cp = 0
            baf_priority = 0
            prof_baf = np.array([])
        
        #if vaf_queue.qsize() > 0:
        if len(vaf_queue) > 0:
            vaf_priority, vaf_clasp_tree_idx = vaf_queue.pop() #vaf_queue.get()
            (vaf_lbound, vaf_ubound), vaf_clasp = vaf_tree[vaf_clasp_tree_idx]
            vaf_cp = vaf_lbound + my_split(vaf_clasp, vaf_clasp.profile, validation=validation, threshold=threshold)
            vaf_profile[vaf_lbound:vaf_ubound - window_size + 1] = np.max([vaf_profile[vaf_lbound:vaf_ubound - window_size + 1], vaf_clasp.profile], axis=0)
            prof_vaf = vaf_clasp.profile
            
        else:
            vaf_cp = 0
            vaf_priority = 0
            prof_vaf = np.array([])
            
        if len(maf_queue) > 0:
            maf_priority, maf_clasp_tree_idx = maf_queue.pop() #vaf_queue.get()
            (maf_lbound, maf_ubound), maf_clasp = maf_tree[maf_clasp_tree_idx]
            maf_cp = maf_lbound + my_split(maf_clasp, maf_clasp.profile, validation=validation, threshold=threshold)
            maf_profile[maf_lbound:vaf_ubound - window_size + 1] = np.max([maf_profile[maf_lbound:maf_ubound - window_size + 1], maf_clasp.profile], axis=0)
            prof_maf = maf_clasp.profile
            
        else:
            maf_cp = 0
            maf_priority = 0
            prof_maf = np.array([])
            
        
        all_bound = [(dr_lbound, dr_ubound), 
                     (baf_lbound, baf_ubound), 
                     (vaf_lbound, vaf_ubound),
                     (maf_lbound, maf_ubound)]
        
        all_cp = [dr_cp, baf_cp, vaf_cp, maf_cp]
        all_score = [dr_priority, baf_priority, vaf_priority, maf_priority]
        
        all_profile = [prof_dr, prof_baf, prof_vaf, prof_maf]
        all_profile = np.array([prof for prof in all_profile if prof.size != 0])
        
                    
        if mode == 'max': # max score
            idx = np.argmax(np.abs(all_score))
            keep_cp = all_cp[idx]
            
            
        elif mode == 'mult':
            new_score = np.prod(all_profile, axis = 0) #prof_dr * prof_baf * prof_vaf
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score)
            keep_cp = dr_lbound + np.argmax(new_score)
            
            
        elif mode == 'sum':
            new_score = np.sum(all_profile, axis = 0)
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score/3)
            keep_cp = dr_lbound + np.argmax(new_score/3)
            
        
        if mode != 'max':
            all_val = [] 
            if prof_vaf.size > 0:
                val_vaf = significance_test(vaf_clasp, test_cp, threshold)
                all_val.append(all_val)
            
            if prof_baf.size > 0:
                val_baf = significance_test(baf_clasp, test_cp, threshold)
                all_val.append(val_baf)
            
            if prof_dr.size > 0:
                val_dr = significance_test(dr_clasp, test_cp, threshold)
                all_val.append(val_dr)
            
            if prof_maf.size > 0:
                val_maf = significance_test(maf_clasp, test_cp, threshold)
                all_val.append(val_maf)
            
            if (True in all_val):
                CP.append(keep_cp)
        
        else:
            CP.append(keep_cp)
            
        
        
        dr_scores.append(-dr_priority)
        baf_scores.append(-baf_priority)
        vaf_scores.append(-vaf_priority)
        maf_scores.append(-maf_priority)
        
        # define new rrange
        lrange, rrange = (dr_lbound, keep_cp), (keep_cp, dr_ubound) 
        for prange in (lrange, rrange):
            #print('\n Search in bound:', prange)
            
            low = prange[0]
            high = prange[1]

            #print('\n DR')
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
            #print('\n BAF')
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

            #print('\n VAF')
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
            
            maf_tree, maf_queue, maf_scores_tmp = local_segmentation(low, 
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
                                                                    maf, 
                                                                    validation, 
                                                                    threshold,
                                                                    maf_tree, 
                                                                    maf_queue,
                                                                    n_timepoints)
            
            #print(f'scores = {dr_scores_tmp}, {baf_scores_tmp}, {vaf_scores_tmp}')
            if dr_scores_tmp != None:
                dr_scores.append(dr_scores_tmp)
                
            if baf_scores_tmp != None:
                baf_scores.append(baf_scores_tmp)
                
            if vaf_scores_tmp != None:
                vaf_scores.append(vaf_scores_tmp)
                
            if maf_scores_tmp != None:
                maf_scores.append(maf_scores_tmp)
    
    return CP
    