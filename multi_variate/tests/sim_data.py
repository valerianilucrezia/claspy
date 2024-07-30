import numpy as np

"""Function to merge simulated distributions of stable and mutant allele frequencies and determine change points. Accessed by sim_data function
Input:
    - baseline: 1darray of randomly generated numbers between x and y of length size_segment representing the frequency of non-mutated bases
    - mutant: 1darray of randomly generated numbers between x and y of length size_segment representing the frequency of mutated bases
    - size_segment: int representing the size of each region being simulated
    - size_ts: Int representing the size of the time series. Should be a multiple of size_segment
    - prob: Tuple representing the probability of a baseline frequency of a mutant frequency
    - seed: Int to set random number generator for reproducibility
Output:
    - profile: concatenated profile of baseline and mutant frequencies
    - cp: list of known change points from simulated profile
"""
def _random_profile_and_cp(baseline:np.ndarray, mutant:np.ndarray, size_segment:int, size_ts:int, prob:tuple, seed:int):
    rng = np.random.Generator(np.random.PCG64(seed=seed))
    expression_sets = rng.choice((baseline, mutant), size=int(size_ts/size_segment), p=prob)

    cp = set()
    for i in range(1, len(expression_sets)-1):
        if np.array_equiv(expression_sets[i], expression_sets[i-1]) is False:
            cp.add(i*size_segment)

    profile = np.concatenate(expression_sets, axis=None)
    return profile, sorted(cp)


"""Function to simulate distributions of stable allele frequencies and mutated allele frequencies and return a known profile time series and change points
Input:
    - base_freq_range: Tuple representing the lower and upper range to sample from for constructing a baseline frequency
    - mutant_freq_range: Tuple representing the lower and upper range to sample from for constructing a mutant frequency
    - size_segment: Int representing the size of each region being simulated. Default 20
    - size_ts: Int representing the size of the time series. Should be a multiple of size_segment. Default 1000
    - prob: Tuple representing the probability of a baseline frequency of a mutant frequency. Default (0.9, 0.1)
    - seed: Int to set random number generator for reproducibility
Output:
    - profile: 2Darray of time series points and profile of baseline and mutant frequencies
    - cp: list of known change points from simulated profile
"""
def sim_data(base_freq_range:tuple, mutant_freq_range:tuple, size_segment:int=20, size_ts:int=1000, prob:tuple=(0.9, 0.1), seed:int=0):
    np.random.seed(seed)
    baseline = np.random.uniform(base_freq_range[0], base_freq_range[1], size=size_segment)
    mutant = np.random.uniform(mutant_freq_range[0], mutant_freq_range[1], size=size_segment)
    profile, cp = _random_profile_and_cp(baseline, mutant, size_segment=size_segment, size_ts=size_ts, prob=prob, seed=seed)
    ts = np.array((range(len(profile)),profile))
    return ts, cp
