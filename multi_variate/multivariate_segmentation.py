from claspy.utils import check_input_time_series, check_excl_radius
from claspy.clasp import ClaSPEnsemble
from claspy.validation import map_validation_tests
from claspy.validation import significance_test
from claspy.segmentation import BinaryClaSPSegmentation
import numpy as np

"""
Class: multivariateClaSP
Inherits from BinaryClaSPSegmentation from claspy package
"""
class MultivariateClaSPSegmentation(BinaryClaSPSegmentation):
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
                
        self.cp = my_split(self.clasp, self.clasp.profile, validation=self.validation, threshold=self.threshold)
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


"""Check that a change point is valid"""
def cp_is_valid(candidate, 
                change_points, 
                n_timepoints, 
                min_seg_size):
    
    for change_point in [0] + change_points + [n_timepoints]:
        left_begin = max(0, change_point - min_seg_size)
        right_end = min(n_timepoints, change_point + min_seg_size)
        
        if candidate in range(left_begin, right_end): 
            return False

    return True


"""This function can return a value of None if the significance test for the change points is not significant, or return the change point
    Successive functions that use my_split check that the return value is/not None"""
def my_split(clasp, 
            profile, 
            sparse=True, 
            validation="significance_test", 
            threshold=1e-15):
    
    cp = np.argmax(profile)
   # print(cp)
    # FIXME what if the validation test is None? And what is sparse?
    if validation is not None:
        validation_test = map_validation_tests(validation)
        # validation test evaluates to True or False, so if false the cp is not significant
        if not validation_test(clasp, cp, threshold):
            return None

    if sparse is True:
        return cp


def take_first_cp(multivariate_clasp_objects: dict, 
                  mode):
    
    profiles = [i.profiles for i in multivariate_clasp_objects.values()]
    cps = [np.argmax(i) for i in profiles]
    scores = [max(i) for i in profiles]

    if mode == 'max':
        idx = np.argmax(scores)
        most_common_cp = cps[idx]
            
    elif mode == 'mult':
        most_common_cp = np.argmax(np.prod(np.array(profiles), axis=0))
        
    elif mode == 'sum':
        most_common_cp = np.argmax(np.sum(np.array(profiles), axis=0)/len(multivariate_clasp_objects.values()))
        
    return most_common_cp


# What is validation test referring to?
def validate_first_cp(multivariate_clasp_objects: dict,
                      cp):
    # validate clasp object for each one given as input
    for ts_obj in multivariate_clasp_objects.values():
        validation_test = map_validation_tests(ts_obj.validation)
        ts_obj.val = validation_test(ts_obj.clasp, cp, ts_obj.threshold)
    # append clasp object variables
    if any([i.val for i in multivariate_clasp_objects.values()]):
        for ts_obj in multivariate_clasp_objects.values():
            ts_obj.tree.append((ts_obj.range, ts_obj.clasp))
            ts_obj.queue.append((-ts_obj.profile[cp], len(ts_obj.tree) - 1))
    
    
# TODO: Proofread this one
def find_cp_iterative(multivariate_clasp_objects, mode):
    CP = []

    # n_segments should be the same for each from what I can tell, so just get one
    n_segments = multivariate_clasp_objects.values()[0].n_segments
    
    for _ in range(n_segments - 1):
        
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
                mysplit_ret = my_split(ts_obj.clasp, ts_obj.clasp.profile, validation=ts_obj.validation, threshold=ts_obj.threshold)
                if mysplit_ret is None:
                    mysplit_ret = 0
                ts_obj.cp = ts_obj.lbound + mysplit_ret
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
    
    return CP
