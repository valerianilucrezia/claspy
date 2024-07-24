# load in text/pkl file of true and pred bps from first iteration of code
from main import main
import pickle as pkl
def load_standards(file):
    standard = pkl.load(file)
    # assume dict
    return (standard["true_bps"], standard["pred_bps"])

# load in original data
# def load data...

# run original data in current main file
# results = main(original_data)

# check if results have changed. If yes, test fails
truth = load_standards
# TODO: can either make main an object or make the entire function of main for regression testing
# making object sounds more reasonable
results = main
# assert truth.true_bps == results.true_bps
# assert truth.pred_bps == results.pred_bps