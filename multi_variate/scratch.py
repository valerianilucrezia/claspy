from tests.sim_data import sim_data, sim_data_baf
from multivariate_clasp import MultivariateClaSP
import numpy as np
import pandas as pd
import os

def flatten(xss):
    return [x for xs in xss for x in xs]

def main():
    baf_profile, baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=400, size_ts=10000, prob_baseline=0.95, seed=97)
    vaf_profile, vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.25, 1), seed=93, size_ts=10000, size_segment=400, prob_baseline=0.95)
    dr_profile, dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=74, size_ts=10000, size_segment=400, prob_baseline=0.95)

    # cna_id field is for checking accuracy of change points. Since change points are simulated and thus known already, cna_id here is array of zeros
    test_df = pd.DataFrame({"vaf":vaf_profile, "median_baf":baf_profile, "median_dr":dr_profile, "pos":list(range(0, len(vaf_profile))), "cna_id":np.zeros(len(vaf_profile))})
    
    # make temporary directory for testing
    tmp = "/Users/brandonhastings/code_projects/claspy/multi_variate/scratch_test"
    if os.path.isdir(tmp) is False:
        os.mkdir(tmp)

    
    # make, save, and read csv
    tmp_csv = os.path.join(tmp, "test_csv.csv")
    test_df.to_csv(path_or_buf=tmp_csv, sep=",")
    # csv = pd.read_csv(filepath_or_buffer=tmp_csv, sep=",")
    multiClasp = MultivariateClaSP(tmp_csv, mode='mult', out_dir=os.path.join(tmp, "out"), threshold=1e-15, window_size='suss')
    multiClasp.analyze_time_series()

    clasp_trees = [multiClasp.multivariate_clasp_objects[freq].clasp_tree for freq in multiClasp.frequencies]
    clasp_trees_flat = flatten(clasp_trees)
    clasp = [i[1] for i in clasp_trees_flat]
    for i in clasp:
        print(i.knn.offsets)
        print(i.knn.offsets.shape)

    for i in multiClasp.frequencies:
        print(multiClasp.multivariate_clasp_objects[i].clasp_tree)


    print(clasp_trees_flat)
    print(clasp)

if __name__ == "__main__":
    main()