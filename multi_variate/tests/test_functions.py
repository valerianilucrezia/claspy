from unittest import TestCase
import numpy as np
import multivariate_segmentation as multivariate_segmentation
from clustering import calculate_cluster_score, construct_confusion_matrix
from sim_data import sim_data, sim_data_dr


# simulate ground truths
vaf_profile, vaf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.9, 0.95))
baf_profile, baf_cp = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.9, 0.95))
dr_profile, dr_cp = sim_data(base_freq_range=(0.8, 1.1), mutant_freq_range=(1.5, 2.1))


class TestMultivariateClaSPSegementation(TestCase):
    def setUp(self) -> None:
        pass

class TestGetFirstCP(TestCase):
    pass

class TestCpIsValid(TestCase):
    pass

class TestMySplit(TestCase):
    pass

class TestLocalSegmentation(TestCase):
    pass

class TestTakeFirstCP(TestCase):
    pass

class TestValidateFirstCP(TestCase):
    pass

class TestFindCPIterative(TestCase):
    pass

class TestMAFFindCPIterative(TestCase):
    pass

class TestCalculateClusterScore(TestCase):
    baseline = np.random.uniform(0.1, 0.2, 20)
    mutant = np.random.uniform(0.5, 0.6, 20)
    profile = np.concatenate((baseline, mutant, baseline, mutant, baseline), axis=None)
    ts = np.array((range(100),profile))
    cp = range(21, 100, 20)
    labels, scores = calculate_cluster_score(ts, cp)
    print(cp)
    print(scores)