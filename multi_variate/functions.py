from queue import PriorityQueue
from claspy.utils import check_input_time_series, check_excl_radius
from claspy.clasp import ClaSPEnsemble
from claspy.validation import map_validation_tests
from claspy.validation import significance_test_2
from claspy.validation import significance_test
from claspy.segmentation import BinaryClaSPSegmentation

import numpy as np

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Class: multivariateClaSP
Inherits from BinaryClaSPSegmentation from claspy package

"""
class multivariateClaSP(BinaryClaSPSegmentation):
    def __init__(self, time_series, n_segments="learn", n_estimators=10, window_size="suss", k_neighbours=3, distance="znormed_euclidean_distance", score="roc_auc", early_stopping=True, validation="significance_test", threshold=1e-15, excl_radius=5, n_jobs=-1, random_state=2357):
        super().__init__(n_segments, n_estimators, window_size, k_neighbours, distance, score, early_stopping, validation, threshold, excl_radius, n_jobs, random_state)
        
        self.time_series = time_series
        self.min_seg_size = self.window_size * self.excl_radius

        # self.min_seg_size = self.window_size * self.excl_radius
    def get_first_cp(self):
        
        check_input_time_series(self.time_series)
        check_excl_radius(self.k_neighbours, self.excl_radius)

        n_timepoints = self.time_series.shape[0]
        min_seg_size = self.window_size * self.excl_radius
        
        
        self.queue = [] #PriorityQueue()
        self.clasp_tree = []
        
        if n_segments == "learn":
            n_segments = self.time_series.shape[0] // min_seg_size
            
        self.prange = 0, self.time_series.shape[0]    
        self.clasp = ClaSPEnsemble(n_estimators=self.n_estimators,
                        window_size=self.window_size,
                        k_neighbours=self.k_neighbours,
                        distance=self.distance,
                        score=self.score,
                        early_stopping=self.early_stopping,
                        excl_radius=self.excl_radius,
                        n_jobs=self.n_jobs,
                        random_state=self.random_state).fit(self.time_series, 
                                                    validation=self.validation, 
                                                    threshold=self.threshold)
                
        self.cp = my_split(self.clasp, self.clasp.profile, validation=self.validation, threshold=self.threshold) #clasp.split(validation=validation, threshold=threshold)
        # print('first cp',cp)
        # # check if it is valid with significant test
        
        self.profile = self.clasp.profile
        return(self.clasp, self.profile, self.cp, self.prange, self.clasp_tree, self.queue)
    
    def local_segmentation(self, lbound, ubound, change_points):

        if ubound - lbound < 2 * self.min_seg_size: 
            return self.clasp_tree, self.queue, 0

        self.clasp = ClaSPEnsemble(
            n_estimators=self.n_estimators,
            window_size=self.window_size,
            k_neighbours=self.k_neighbours,
            distance=self.distance,
            score=self.score,
            early_stopping=self.early_stopping,
            excl_radius=self.excl_radius,
            n_jobs=self.n_jobs,
            random_state=self.random_state
        ).fit(self.time_series[lbound:ubound], validation=self.validation, threshold=self.threshold)
        
        cp = my_split(self.clasp, self.clasp.profile, validation=self.validation, threshold=self.threshold)
        # FIXME: cp will never be None
        if cp is None: 
            return self.clasp_tree, self.queue, 0
        
        # FIXME: score is originally a string argument and here (and below) it is reclassified as int. Change name?
        self.score = self.clasp.profile[cp]
        #print('cp, score=', cp, score)
        
        if not cp_is_valid(lbound + cp, change_points, self.n_timepoints, self.min_seg_size):  #candidate, change_points, n_timepoints, min_seg_size
            #print('cp is not valid')
            #return clasp_tree, queue, 0
            self.score = 0

        self.clasp_tree.append(((lbound, ubound), self.clasp))
        #queue.put((-score, len(clasp_tree) - 1))    
        self.queue.append((-self.score, len(self.clasp_tree) - 1))
        return self.clasp_tree, self.queue, self.score



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

# FIXME: this function either returns an int or a tuple of vectors. Sussecive functions check for varying return types that can then fail
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
            # TODO: Why not None? local_segmentation checks if None, but find_cp_iterative uses it as int

            return 0 #None

    if sparse is True:
        return cp
    # FIXME: why does this return a variable that is not created?
    return time_series[:cp], time_series[cp:]

# ----- Moved to MultivariateClaSP class ------

# def local_segmentation(lbound, 
#                        ubound, 
#                        change_points, 
#                        min_seg_size,
#                        n_estimators,
#                        window_size,
#                        k_neighbours,
#                        distance, 
#                        score,
#                        early_stopping,
#                        excl_radius,
#                        n_jobs,
#                        random_states,
#                        time_series,
#                        validation,
#                        threshold,
#                        clasp_tree,
#                        queue, 
#                        n_timepoints):

#     if ubound - lbound < 2 * min_seg_size: 
#         return clasp_tree, queue, 0

#     clasp = ClaSPEnsemble(
#         n_estimators=n_estimators,
#         window_size=window_size,
#         k_neighbours=k_neighbours,
#         distance=distance,
#         score=score,
#         early_stopping=early_stopping,
#         excl_radius=excl_radius,
#         n_jobs=n_jobs,
#         random_state=random_states
#     ).fit(time_series[lbound:ubound], validation=validation, threshold=threshold)
    
#     cp = my_split(clasp, clasp.profile, validation=validation, threshold=threshold)
    
#     if cp is None: 
#         return clasp_tree, queue, 0
    
#     score = clasp.profile[cp]
#     #print('cp, score=', cp, score)
    
#     if not cp_is_valid(lbound + cp, change_points, n_timepoints, min_seg_size):  #candidate, change_points, n_timepoints, min_seg_size
#         #print('cp is not valid')
#         #return clasp_tree, queue, 0
#         score = 0

#     clasp_tree.append(((lbound, ubound), clasp))
#     #queue.put((-score, len(clasp_tree) - 1))    
#     queue.append((-score, len(clasp_tree) - 1))
#     return clasp_tree, queue, score

# FIXME: Needs refactored to be method of multivariateClaSP to add a variable .first_cp
def take_first_cp(dr_profile, 
                  vaf_profile, 
                  baf_profile, 
                  mode, 
                  maf_profile = None):    
    
    if maf_profile is not None:
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
    
# TODO: should a MAF object be created in the try loop in main?
# What is validation test referring to?
def validate_first_cp(multivariate_clasp_objects, cp):
    
    # validate clasp object for each one given as input
    for ts_obj in multivariate_clasp_objects:
        validation_test = map_validation_tests(ts_obj.validation)
        ts_obj.val = validation_test(ts_obj.clasp, cp, ts_obj.threshold)
    # append clasp object variables
    if any([i.val for i in multivariate_clasp_objects]):
        for ts_obj in multivariate_clasp_objects:
            ts_obj.tree.append((ts_obj.range, ts_obj.clasp))
            ts_obj.queue.append((-ts_obj.profile[cp], len(ts_obj.tree) - 1))
    
    
# TODO: appears stable, but could be a method of MultivariateClaSP
def find_cp_iterative(multivariate_clasp_objects, mode):
    CP = []
    dr_scores = []
    baf_scores = []
    vaf_scores = []

    # n_segments should be the same for each from what I can tell, so just get one
    n_segments = multivariate_clasp_objects.values()[0].n_segments
    
    for idx in range(n_segments - 1):
        
        if all([i.queue for i in multivariate_clasp_objects.values()], 0):
            print('Stop because queue is empty')
            break
        # TODO: Make sure it is okay that the class variables are modified, pretty sure it's fine
        # TODO: is ts_obj.profile, .clasp.profile, and now .prof all the same thing?
        for ts_obj in multivariate_clasp_objects.values():
            if len(ts_obj.queue) > 0:
                ts_obj.priority, ts_obj.clasp_tree_idx = ts_obj.queue.pop() #dr_queue.get()
                (ts_obj.lbound, ts_obj.ubound), ts_obj.clasp = ts_obj.tree[ts_obj.clasp_tree_idx]
                # FIXME: Here if mysplit returns None an error would occur. Possible solution below. Alternative: try except
                # mysplit_ret = my_split(ts_obj.clasp, ts_obj.clasp.profile, validation=ts_obj.validation, threshold=ts_obj.threshold)
                # if mysplit_ret is None:
                #     mysplit_ret = 0
                # ts_obj.cp = ts_obj.lbound + mysplit_ret
                ts_obj.cp = ts_obj.lbound + my_split(ts_obj.clasp, ts_obj.clasp.profile, validation=ts_obj.validation, threshold=ts_obj.threshold)
                ts_obj.profile[ts_obj.lbound:ts_obj.ubound - ts_obj.window_size + 1] = np.max([ts_obj.profile[ts_obj.lbound:ts_obj.ubound - ts_obj.window_size + 1], ts_obj.clasp.profile], axis=0)
                ts_obj.prof = ts_obj.clasp.profile
            
            else:
                ts_obj.cp = 0
                ts_obj.priority = 0
                ts_obj.prof = np.array([])
            
        all_bound = [(i.lbound, i.ubound) for i in multivariate_clasp_objects.values()]
        
        all_cp = [i.cp for i in multivariate_clasp_objects.values()]
        all_score = [i.priority for i in multivariate_clasp_objects.values()]
        
        all_profile = np.array([i.prof for i in multivariate_clasp_objects.values() if i.prof.size != 0])
        
        # TODO: What do these do
        if mode == 'max': # max score
            keep_cp = all_cp[np.argmax(np.abs(all_score))]
            
        elif mode == 'mult':
            new_score = np.prod(all_profile, axis = 0) #prof_dr * prof_baf * prof_vaf
            # TODO: Need test case
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score)
            keep_cp = multivariate_clasp_objects["dr"].lbound + np.argmax(new_score)
            
            
        elif mode == 'sum':
            new_score = np.sum(all_profile, axis = 0)
            new_score[new_score == np.inf] = -10000
            test_cp = np.argmax(new_score/3)
            keep_cp = multivariate_clasp_objects["dr"].lbound + np.argmax(new_score/3)
            
        
        if mode != 'max':
            all_val = [] 
            for ts_obj in multivariate_clasp_objects.values():
                if ts_obj.prof.size > 0:
                    val = significance_test(ts_obj.clasp, test_cp, ts_obj.threshold)
                    all_val.append(val)
                if (True in all_val):
                    CP.append(keep_cp)
        
        else:
            CP.append(keep_cp)
            
        # FIXME: If we are only appending CP, nothing below seems to matter
        # FIXME: Can append the priority variable of the ts object instead
        # dr_scores.append(-dr_priority)
        # baf_scores.append(-baf_priority)
        # vaf_scores.append(-vaf_priority)
        
        # define new rrange
        lrange, rrange = (multivariate_clasp_objects["dr"].lbound, keep_cp), (keep_cp, multivariate_clasp_objects["dr"].ubound) 
        for prange in (lrange, rrange):
            #print('\n Search in bound:', prange)
            
            low = prange[0]
            high = prange[1]
            for ts_obj in multivariate_clasp_objects.values():
                ts_obj.local_segmentation(low, high, CP)

            # FIXME: refactored, but why is this here if it is not returning the *_scores lists in the end?
            # print(f'scores = {dr_scores_tmp}, {baf_scores_tmp}, {vaf_scores_tmp}')
            # if dr_scores_tmp != None:
            #     dr_scores.append(dr_scores_tmp)
                
            # if baf_scores_tmp != None:
            #     baf_scores.append(baf_scores_tmp)
                
            # if vaf_scores_tmp != None:
            #     vaf_scores.append(vaf_scores_tmp)
            # FIXME: Should a score object be created as above and appended to?
    
    return CP

# TODO: I think I have made this function unnecesary
# def MAF_find_cp_iterative(dr_clasp, dr_tree, dr_queue, dr_profile,
#                       baf_clasp, baf_tree, baf_queue, baf_profile,
#                       vaf_clasp, vaf_tree, vaf_queue, vaf_profile,
#                       n_segments, validation, threshold, window_size, min_seg_size,
#                       n_estimators, k_neighbours, distance, scored, early_stopping, 
#                       excl_radius, n_jobs, random_state, n_timepoints, 
#                       dr, baf, vaf, 
#                       mode,
#                       maf_clasp, maf_tree, maf_queue, maf_profile, maf):
#     CP = []
#     dr_scores = []
#     baf_scores = []
#     vaf_scores = []
#     maf_scores = []
    
#     for idx in range(n_segments - 1):
#         if len(dr_queue) == 0 and len(baf_queue) == 0 and len(vaf_queue) == 0 and len(maf_queue) == 0:
#             print('Stop because queue is empty')
#             break
        
#         if len(dr_queue) > 0:
#             dr_priority, dr_clasp_tree_idx = dr_queue.pop() #dr_queue.get()
#             (dr_lbound, dr_ubound), dr_clasp = dr_tree[dr_clasp_tree_idx]
#             # FIXME: Here if mysplit returns None an error would occur
#             dr_cp = dr_lbound + my_split(dr_clasp, dr_clasp.profile, validation=validation, threshold=threshold)
#             dr_profile[dr_lbound:dr_ubound - window_size + 1] = np.max([dr_profile[dr_lbound:dr_ubound - window_size + 1], dr_clasp.profile], axis=0)
#             prof_dr = dr_clasp.profile
        
#         else:
#             dr_cp = 0
#             dr_priority = 0
#             prof_dr = np.array([])
        

#         #if baf_queue.qsize() > 0:
#         if len(baf_queue) > 0: 
#             baf_priority, baf_clasp_tree_idx = baf_queue.pop() #baf_queue.get()
#             (baf_lbound, baf_ubound), baf_clasp = baf_tree[baf_clasp_tree_idx]
#             baf_cp = baf_lbound + my_split(baf_clasp, baf_clasp.profile, validation=validation, threshold=threshold)
#             baf_profile[baf_lbound:baf_ubound - window_size + 1] = np.max([baf_profile[baf_lbound:baf_ubound - window_size + 1], baf_clasp.profile], axis=0)
#             prof_baf = baf_clasp.profile
        
#         else:
#             baf_cp = 0
#             baf_priority = 0
#             prof_baf = np.array([])
        
#         #if vaf_queue.qsize() > 0:
#         if len(vaf_queue) > 0:
#             vaf_priority, vaf_clasp_tree_idx = vaf_queue.pop() #vaf_queue.get()
#             (vaf_lbound, vaf_ubound), vaf_clasp = vaf_tree[vaf_clasp_tree_idx]
#             vaf_cp = vaf_lbound + my_split(vaf_clasp, vaf_clasp.profile, validation=validation, threshold=threshold)
#             vaf_profile[vaf_lbound:vaf_ubound - window_size + 1] = np.max([vaf_profile[vaf_lbound:vaf_ubound - window_size + 1], vaf_clasp.profile], axis=0)
#             prof_vaf = vaf_clasp.profile
            
#         else:
#             vaf_cp = 0
#             vaf_priority = 0
#             prof_vaf = np.array([])
            
#         if len(maf_queue) > 0:
#             maf_priority, maf_clasp_tree_idx = maf_queue.pop() #vaf_queue.get()
#             (maf_lbound, maf_ubound), maf_clasp = maf_tree[maf_clasp_tree_idx]
#             maf_cp = maf_lbound + my_split(maf_clasp, maf_clasp.profile, validation=validation, threshold=threshold)
#             maf_profile[maf_lbound:vaf_ubound - window_size + 1] = np.max([maf_profile[maf_lbound:maf_ubound - window_size + 1], maf_clasp.profile], axis=0)
#             prof_maf = maf_clasp.profile
            
#         else:
#             maf_cp = 0
#             maf_priority = 0
#             prof_maf = np.array([])
            
        
#         all_bound = [(dr_lbound, dr_ubound), 
#                      (baf_lbound, baf_ubound), 
#                      (vaf_lbound, vaf_ubound),
#                      (maf_lbound, maf_ubound)]
        
#         all_cp = [dr_cp, baf_cp, vaf_cp, maf_cp]
#         all_score = [dr_priority, baf_priority, vaf_priority, maf_priority]
        
#         all_profile = [prof_dr, prof_baf, prof_vaf, prof_maf]
#         all_profile = np.array([prof for prof in all_profile if prof.size != 0])
        
                    
#         if mode == 'max': # max score
#             idx = np.argmax(np.abs(all_score))
#             keep_cp = all_cp[idx]
            
            
#         elif mode == 'mult':
#             new_score = np.prod(all_profile, axis = 0) #prof_dr * prof_baf * prof_vaf
#             new_score[new_score == np.inf] = -10000
#             test_cp = np.argmax(new_score)
#             keep_cp = dr_lbound + np.argmax(new_score)
            
            
#         elif mode == 'sum':
#             new_score = np.sum(all_profile, axis = 0)
#             new_score[new_score == np.inf] = -10000
#             test_cp = np.argmax(new_score/3)
#             keep_cp = dr_lbound + np.argmax(new_score/3)
            
        
#         if mode != 'max':
#             all_val = [] 
#             if prof_vaf.size > 0:
#                 val_vaf = significance_test(vaf_clasp, test_cp, threshold)
#                 all_val.append(all_val)
            
#             if prof_baf.size > 0:
#                 val_baf = significance_test(baf_clasp, test_cp, threshold)
#                 all_val.append(val_baf)
            
#             if prof_dr.size > 0:
#                 val_dr = significance_test(dr_clasp, test_cp, threshold)
#                 all_val.append(val_dr)
            
#             if prof_maf.size > 0:
#                 val_maf = significance_test(maf_clasp, test_cp, threshold)
#                 all_val.append(val_maf)
            
#             if (True in all_val):
#                 CP.append(keep_cp)
        
#         else:
#             CP.append(keep_cp)
            
        
        
#         dr_scores.append(-dr_priority)
#         baf_scores.append(-baf_priority)
#         vaf_scores.append(-vaf_priority)
#         maf_scores.append(-maf_priority)
        
#         # define new rrange
#         lrange, rrange = (dr_lbound, keep_cp), (keep_cp, dr_ubound) 
#         for prange in (lrange, rrange):
#             #print('\n Search in bound:', prange)
            
#             low = prange[0]
#             high = prange[1]

#             #print('\n DR')
#             dr_tree, dr_queue, dr_scores_tmp = local_segmentation(low, 
#                                                                 high, 
#                                                                 CP, 
#                                                                 min_seg_size,
#                                                                 n_estimators,
#                                                                 window_size,
#                                                                 k_neighbours,
#                                                                 distance, 
#                                                                 scored,
#                                                                 early_stopping,
#                                                                 excl_radius,
#                                                                 n_jobs,
#                                                                 random_state,
#                                                                 dr, 
#                                                                 validation, 
#                                                                 threshold,
#                                                                 dr_tree, 
#                                                                 dr_queue,
#                                                                 n_timepoints)
#             #print('\n BAF')
#             baf_tree, baf_queue, baf_scores_tmp = local_segmentation(low, 
#                                                                     high, 
#                                                                     CP, 
#                                                                     min_seg_size,
#                                                                     n_estimators,
#                                                                     window_size,
#                                                                     k_neighbours,
#                                                                     distance, 
#                                                                     scored,
#                                                                     early_stopping,
#                                                                     excl_radius,
#                                                                     n_jobs,
#                                                                     random_state,
#                                                                     baf, 
#                                                                     validation, 
#                                                                     threshold,
#                                                                     baf_tree, 
#                                                                     baf_queue,
#                                                                     n_timepoints)

#             #print('\n VAF')
#             vaf_tree, vaf_queue, vaf_scores_tmp = local_segmentation(low, 
#                                                                     high, 
#                                                                     CP, 
#                                                                     min_seg_size,
#                                                                     n_estimators,
#                                                                     window_size,
#                                                                     k_neighbours,
#                                                                     distance, 
#                                                                     scored,
#                                                                     early_stopping,
#                                                                     excl_radius,
#                                                                     n_jobs,
#                                                                     random_state,
#                                                                     vaf, 
#                                                                     validation, 
#                                                                     threshold,
#                                                                     vaf_tree, 
#                                                                     vaf_queue,
#                                                                     n_timepoints)
            
#             maf_tree, maf_queue, maf_scores_tmp = local_segmentation(low, 
#                                                                     high, 
#                                                                     CP, 
#                                                                     min_seg_size,
#                                                                     n_estimators,
#                                                                     window_size,
#                                                                     k_neighbours,
#                                                                     distance, 
#                                                                     scored,
#                                                                     early_stopping,
#                                                                     excl_radius,
#                                                                     n_jobs,
#                                                                     random_state,
#                                                                     maf, 
#                                                                     validation, 
#                                                                     threshold,
#                                                                     maf_tree, 
#                                                                     maf_queue,
#                                                                     n_timepoints)
            
#             #print(f'scores = {dr_scores_tmp}, {baf_scores_tmp}, {vaf_scores_tmp}')
#             if dr_scores_tmp != None:
#                 dr_scores.append(dr_scores_tmp)
                
#             if baf_scores_tmp != None:
#                 baf_scores.append(baf_scores_tmp)
                
#             if vaf_scores_tmp != None:
#                 vaf_scores.append(vaf_scores_tmp)
                
#             if maf_scores_tmp != None:
#                 maf_scores.append(maf_scores_tmp)
    
#     return CP
    