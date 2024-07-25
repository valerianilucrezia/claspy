from unittest import TestCase
from main import MultivariateClaSP
import pickle as pkl

"""check if results (cp) have changed from original MVP run. If changed, test fails"""
class TestRegression(TestCase):
    # settings (change later)
    mode = None
    genes_file = None
    standards_file = None

    truth = pkl.load(standards_file)

    # FIXME fix this in main so it doesn't save anything out - in progress
    results = MultivariateClaSP(genes_file, mode=mode, out_dir=None)
    results.analyze_time_series()

    assert truth == results.CP
