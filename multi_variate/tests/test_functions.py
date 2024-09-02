import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from unittest import TestCase, main
import numpy as np
import pandas as pd
import tempfile
from multivariate_clasp import MultivariateClaSP
from clustering import cluster_segments, calculate_cluster_scores, construct_confusion_matrix
from sim_data import sim_data, sim_data_baf
from data_import import get_data_csv, get_data_tsv


# simulate ground truths
# baf_profile, baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=2000, size_ts=10000, prob_baseline=0.95, seed=11)
# vaf_profile, vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.1, 0.105), seed=4, size_ts=10000, size_segment=2000, prob_baseline=0.95)
# dr_profile, dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=1, size_ts=10000, size_segment=2000, prob_baseline=0.95)


class TestDataImport(TestCase):
 
    def setUp(self):

        # simulate ground truths
        self.baf_profile, self.baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=2000, size_ts=10000, prob_baseline=0.95, seed=11)
        self.vaf_profile, self.vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.1, 0.105), seed=4, size_ts=10000, size_segment=2000, prob_baseline=0.95)
        self.dr_profile, self.dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=1, size_ts=10000, size_segment=2000, prob_baseline=0.95)

        # convert ground truths to dataframe
        self.test_df = pd.DataFrame({"vaf":self.vaf_profile, "median_baf":self.baf_profile, "median_dr":self.dr_profile, "pos":list(range(0, len(self.vaf_profile))), "cna_id":np.zeros(len(self.vaf_profile))}).round(8)

        # make temporary directory for testing
        self.tmp = tempfile.TemporaryDirectory()
        
        # make csv path
        self.tmp_csv = os.path.join(self.tmp.name, "test_csv")
        # save csv
        self.test_df.to_csv(path_or_buf=self.tmp_csv, sep=",", float_format='%.8f')
        # get needed data into dict
        self.res_dict = get_data_csv(data=self.tmp_csv, frequencies=["vaf", "median_baf", "median_dr"])

        # # make, save, and read tsv
        # tmp_tsv = os.path.join(self.tmp.name, "test_tsv")
        # self.test_df.to_csv(path_or_buf=tmp_tsv, sep="\t")
        # self.tsv = pd.read_csv(filepath_or_buffer=tmp_tsv, sep=",")

    def test_baf_import(self):
        self.assertTrue(np.allclose(a=self.res_dict["median_baf"], b=self.baf_profile))
        self.assertTrue(type(self.res_dict["median_baf"])==type(self.baf_profile))
        self.assertTrue(len(self.res_dict["median_baf"])==len(self.baf_profile))
    
    def test_dr_import(self):
        self.assertTrue(np.allclose(a=self.res_dict["median_dr"], b=self.dr_profile))
        self.assertTrue(type(self.res_dict["median_dr"])==type(self.dr_profile))
        self.assertTrue(len(self.res_dict["median_dr"])==len(self.dr_profile))
    
    def test_vaf_import(self):
        self.assertTrue(np.allclose(a=self.res_dict["vaf"], b=self.vaf_profile))
        self.assertTrue(type(self.res_dict["vaf"])==type(self.vaf_profile))
        self.assertTrue(len(self.res_dict["vaf"])==len(self.vaf_profile))

    def test_warning_import(self):
        with self.assertWarns(Warning) as cm:
            get_data_csv(data=self.tmp_csv, frequencies=["test"])
        the_warning = cm.warning
        self.assertEqual(type(the_warning), Warning)

    def test_empty_dict_error(self):
        with self.assertRaises(RuntimeError) as cm:
            get_data_csv(data=self.tmp_csv, frequencies=["test"])
        the_exception = cm.exception
        self.assertEqual(type(the_exception), RuntimeError)
        

    def tearDown(self):
        self.tmp.cleanup()

        
class TestMultivariateClaSP_2Breakpoints_2000Segment_1e15Significance(TestCase):
    
    def setUp(self):
        # simulate ground truths
        self.baf_profile, self.baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=2000, size_ts=10000, prob_baseline=0.95, seed=11)
        self.vaf_profile, self.vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.1, 0.105), seed=4, size_ts=10000, size_segment=2000, prob_baseline=0.95)
        self.dr_profile, self.dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=1, size_ts=10000, size_segment=2000, prob_baseline=0.95)

        # cna_id field is for checking accuracy of change points. Since change points are simulated and thus known already, cna_id here is array of zeros
        test_df = pd.DataFrame({"vaf":self.vaf_profile, "median_baf":self.baf_profile, "median_dr":self.dr_profile, "pos":list(range(0, len(self.vaf_profile))), "cna_id":np.zeros(len(self.vaf_profile))})
        
        # make temporary directory for testing
        self.tmp = tempfile.TemporaryDirectory()
        
        # make, save, and read csv
        tmp_csv = os.path.join(self.tmp.name, "test_csv")
        test_df.to_csv(path_or_buf=tmp_csv, sep=",")
        # self.csv = pd.read_csv(filepath_or_buffer=tmp_csv, sep=",")
        self.multiClasp = MultivariateClaSP(tmp_csv, mode='mult', out_dir=os.path.join(self.tmp.name, "out"), threshold=1e-15, window_size='suss')
        self.multiClasp.analyze_time_series()
    
    def test_vaf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points)) == len(self.vaf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points), self.vaf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['vaf'].change_points))} != {len(self.vaf_cp)}")
    
    def test_baf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points)) == len(self.baf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points), self.baf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_baf'].change_points))} != {len(self.baf_cp)}")

    def test_dr_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points)) == len(self.dr_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points), self.dr_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_dr'].change_points))} != {len(self.dr_cp)}")

    def tearDown(self):
        self.tmp.cleanup()



class TestMultivariateClaSP_2Breakpoints_400Segment_1e15Significance(TestCase):
    
    def setUp(self):
        # simulate ground truths
        self.baf_profile, self.baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=400, size_ts=10000, prob_baseline=0.95, seed=99)
        self.vaf_profile, self.vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.1, 0.105), seed=99, size_ts=10000, size_segment=400, prob_baseline=0.95)
        self.dr_profile, self.dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=98, size_ts=10000, size_segment=400, prob_baseline=0.95)

        # cna_id field is for checking accuracy of change points. Since change points are simulated and thus known already, cna_id here is array of zeros
        test_df = pd.DataFrame({"vaf":self.vaf_profile, "median_baf":self.baf_profile, "median_dr":self.dr_profile, "pos":list(range(0, len(self.vaf_profile))), "cna_id":np.zeros(len(self.vaf_profile))})
        
        # make temporary directory for testing
        self.tmp = tempfile.TemporaryDirectory()
        
        # make, save, and read csv
        tmp_csv = os.path.join(self.tmp.name, "test_csv")
        test_df.to_csv(path_or_buf=tmp_csv, sep=",")
        # self.csv = pd.read_csv(filepath_or_buffer=tmp_csv, sep=",")
        self.multiClasp = MultivariateClaSP(tmp_csv, mode='mult', out_dir=os.path.join(self.tmp.name, "out"), threshold=1e-15, window_size='suss')
        self.multiClasp.analyze_time_series()
    
    def test_vaf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points)) == len(self.vaf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points), self.vaf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['vaf'].change_points))} != {len(self.vaf_cp)}")
    
    def test_baf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points)) == len(self.baf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points), self.baf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_baf'].change_points))} != {len(self.baf_cp)}")

    def test_dr_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points)) == len(self.dr_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points), self.dr_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_dr'].change_points))} != {len(self.dr_cp)}")

    def tearDown(self):
        self.tmp.cleanup()


class TestMultivariateClaSP_6Breakpoints_400Segment_1e5Significance(TestCase):
    
    def setUp(self):
        # simulate ground truths
        self.baf_profile, self.baf_cp = sim_data_baf(ploidy=(2,1), purity=1, size_segment=400, size_ts=10000, prob_baseline=0.95, seed=97)
        self.vaf_profile, self.vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.1, 0.105), seed=93, size_ts=10000, size_segment=400, prob_baseline=0.95)
        self.dr_profile, self.dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1), seed=74, size_ts=10000, size_segment=400, prob_baseline=0.95)

        # cna_id field is for checking accuracy of change points. Since change points are simulated and thus known already, cna_id here is array of zeros
        test_df = pd.DataFrame({"vaf":self.vaf_profile, "median_baf":self.baf_profile, "median_dr":self.dr_profile, "pos":list(range(0, len(self.vaf_profile))), "cna_id":np.zeros(len(self.vaf_profile))})
        
        # make temporary directory for testing
        self.tmp = tempfile.TemporaryDirectory()
        
        # make, save, and read csv
        tmp_csv = os.path.join(self.tmp.name, "test_csv")
        test_df.to_csv(path_or_buf=tmp_csv, sep=",")
        # self.csv = pd.read_csv(filepath_or_buffer=tmp_csv, sep=",")
        self.multiClasp = MultivariateClaSP(tmp_csv, mode='mult', out_dir=os.path.join(self.tmp.name, "out"), threshold=1e-5, window_size='suss')
        self.multiClasp.analyze_time_series()
    
    def test_vaf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points)) == len(self.vaf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["vaf"].change_points), self.vaf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['vaf'].change_points))} != {len(self.vaf_cp)}")
    
    def test_baf_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points)) == len(self.baf_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_baf"].change_points), self.baf_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_baf'].change_points))} != {len(self.baf_cp)}")

    def test_dr_cp_detection(self):
        if len(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points)) == len(self.dr_cp):
            abs_vals = np.abs(np.subtract(list(self.multiClasp.multivariate_clasp_objects["median_dr"].change_points), self.dr_cp))
            self.assertTrue(all(i < 50 for i in abs_vals), msg=f"Absolute values of list differences: {abs_vals}")
        else:
            self.fail(msg=f"Assertion Error: {len(list(self.multiClasp.multivariate_clasp_objects['median_dr'].change_points))} != {len(self.dr_cp)}")

    def tearDown(self):
        self.tmp.cleanup()

class TestClustering(TestCase):

    def setUp(self):
        self.profile, self.change_points = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.9, 0.95), size_segment=10, size_ts=100, prob_baseline=0.7)
        # change points for the above are [40, 60, 70, 80, 90], so hand-make the labels to test against
        self.labels = np.concatenate((np.zeros(shape=(40,)),
                                        np.ones(shape=(20,)),
                                        np.zeros(shape=(10,)),
                                        np.ones(shape=(10,)),
                                        np.zeros(shape=(10,)),
                                        np.ones(shape=(10,))), axis=None)

    def test_cluster_segments(self):
        # self.assertListEqual(cluster_segments(list(self.profile), self.change_points), self.labels)
        self.assertTrue(np.array_equal(cluster_segments(self.profile, self.change_points), self.labels))

#     def test_calculate_cluster_scores(self):
#         pass

#     def test_construct_confusion_matrix(self):
#         pass

main()