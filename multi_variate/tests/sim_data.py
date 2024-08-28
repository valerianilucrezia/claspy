import numpy as np


def _random_profile_and_cp(baseline:np.ndarray, mutant:np.ndarray|list, size_segment:int, size_ts:int, prob_baseline:float, seed:int):
    """Function to merge simulated distributions of stable and mutant allele frequencies and determine change points. Accessed by sim_data function.

    Args:
        baseline (np.ndarray): 1darray of randomly generated numbers between x and y of length size_segment representing the frequency of non-mutated bases.

        mutant (np.ndarray | list): 1darray of randomly generated numbers between x and y of length size_segment representing the frequency of mutated bases.

        size_segment (int): The size of each region being simulated.

        size_ts (int): The size of the time series. Should be a multiple of size_segment.

        prob_baseline (float): The probability of a baseline frequency.

        seed (int): Sets the random number generator for reproducibility.

    Returns:
        profile (np.ndarray): 2Darray of time series points and profile of baseline and mutant frequencies.

        cp (np.ndarray): Array of known change points from simulated profile.
    """
    
    rng = np.random.Generator(np.random.PCG64(seed=seed))

    if type(mutant) == list:
        choices = [i for i in mutant]
        choices.append(baseline)

        prob_mutant = (1 - prob_baseline) / len(mutant)
        prob = [prob_mutant for _ in mutant]
        prob.append(prob_baseline)
    
    else:
        choices = (baseline, mutant)
        prob = (prob_baseline, 1 - prob_baseline)

    expression_sets = rng.choice(choices, size=int(size_ts/size_segment), p=prob)

    cp = set()
    for i in range(1, len(expression_sets)):
        if np.array_equiv(expression_sets[i], expression_sets[i-1]) is False:
            cp.add(i*size_segment)

    profile = np.concatenate(expression_sets, axis=None)
    return profile, np.array(sorted(cp))


def sim_data(base_freq_range:tuple, mutant_freq_range:tuple, size_segment:int=20, size_ts:int=1000, prob_baseline:float=0.9, seed:int=0):
    """Function to simulate distributions of stable allele frequencies and mutated allele frequencies and return a known profile time series and change points
    
    Args:
        base_freq_range (tuple): Range representing the lower and upper range to sample from for constructing a baseline frequency.

        mutant_freq_range (tuple): Representing the lower and upper range to sample from for constructing a mutant frequency.

        size_segment (int, optional): Representing the size of each region being simulated. Defaults to 20.
        
        size_ts (int, optional): Representing the size of the time series. Should be a multiple of size_segment. Defaults to 1000.

        prob_baseline (float, optional): The probability of a baseline frequency of a mutant frequency. Defaults to 0.9.

        seed (int, optional): Sets the random number generator for reproducibility. Defaults to 0.
    
    Returns:
        ts (np.ndarray): 2Darray of time series points and profile of baseline and mutant frequencies.
        
        cp (np.ndarray): Array of known change points from simulated profile.
    """
    
    baseline = np.random.uniform(base_freq_range[0], base_freq_range[1], size=size_segment)
    mutant = np.random.uniform(mutant_freq_range[0], mutant_freq_range[1], size=size_segment)
    profile, cp = _random_profile_and_cp(baseline, mutant, size_segment=size_segment, size_ts=size_ts, prob_baseline=prob_baseline, seed=seed)
    ts = np.array(profile)
    return ts, cp


def _calculate_baf(ploidy:tuple=(2,1), purity:float|int=1):
    """Calculate the BAF given a ploidy and purity value.

    Args:
        ploidy (tuple, optional): Ploidy value. Must be tuple of length 2. Defaults to (2,1).

        purity (float | int, optional): Sample purity value between 0 and 1. Defaults to 1.

    Raises:
        ValueError: If purity is not between 0 (exclusive) and 1 (inclusive)

    Returns:
        (tuple): Frequency range to sample a simulated BAF from. Tuple of length 2.
    """
    if not 0 < float(purity) <= 1:
        raise ValueError(f"purity must be between 0 (exclusive) and 1 (inclusive), given {purity}")
    
    frequencies = []
    for n in ploidy:
        baf = (n * purity + (1 - purity)) / ((ploidy[0] + ploidy[1]) * purity + 2 * (1 - purity))
        frequencies.append(baf)
    
    # print(frequencies)
    return tuple(frequencies)


def _construct_baf_segment(frequencies:tuple, size_segment:int, seed:int):
    """Construct a simulated BAF segment of a time series between the given frequencies of length size_segment.

    Args:
        frequencies (tuple): Frequency range to sample a simulated BAF from. Tuple of length 2.

        size_segment (int): Representing the size of each region being simulated.

        seed (int): Sets the random number generator for reproducibility.

    Returns:
        (nd.array): Simulated BAF segment.
    """
    rng = np.random.Generator(np.random.PCG64(seed=seed))
    freq = rng.choice(tuple(frequencies), size=int(size_segment), p=(0.5, 0.5))

    baf_segment = [np.random.uniform(i-0.025, i+0.025) for i in freq]

    return np.array(baf_segment)


def sim_data_baf(ploidy:list|tuple=(2,1), purity:float|int=1, size_segment:int=20, size_ts:int=1000, prob_baseline:float=0.9, seed:int=0):
    """Simulate a BAF time series with mutational stretches equal to length of size segment.

    Args:
        ploidy (tuple, optional): Ploidy value. Must be tuple of length 2. Defaults to (2,1).

        purity (float | int, optional): Sample purity value between 0 and 1. Defaults to 1.

        size_segment (int, optional): Representing the size of each region being simulated. Defaults to 20.
        
        size_ts (int, optional): Representing the size of the time series. Should be a multiple of size_segment. Defaults to 1000.

        prob_baseline (float, optional): The probability of a baseline frequency of a mutant frequency. Defaults to 0.9.

        seed (int, optional): Sets the random number generator for reproducibility. Defaults to 0.

    Returns:
        profile (np.ndarray): 2Darray of time series points and profile of baseline and mutant frequencies.

        cp (np.ndarray): Array of known change points from simulated profile.
    """
    # check this, why 0.5, 0.5 for baseline ploidy?
    # baseline_baf = calculate_baf(ploidy=(0.5, 0.5))
    # shouldn't have to go through calculate baf to construct a baseline
    baseline_baf = (0.5, 0.5)
    # TODO make argument checking better
    if type(ploidy) == list and all(isinstance(i, tuple) for i in ploidy) and type(ploidy[0][0]) == int:
        mutant_bafs = [_calculate_baf(ploidy=i, purity=purity) for i in ploidy]
    elif type(ploidy) == tuple and len(ploidy) == 2 and type(ploidy[0]) == int:
        mutant_bafs = [_calculate_baf(ploidy=ploidy, purity=purity)]

    baseline = _construct_baf_segment(frequencies=baseline_baf, size_segment=size_segment, seed=seed)
    mutants = [_construct_baf_segment(frequencies=i, size_segment=size_segment, seed=seed) for i in mutant_bafs]

    profile, cp = _random_profile_and_cp(baseline=baseline, mutant=mutants, size_segment=size_segment, size_ts=size_ts, prob_baseline=prob_baseline, seed=seed)
    ts = np.array(profile)
    return ts, cp
