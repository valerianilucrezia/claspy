from scipy.stats import ranksums
import numpy as np


def significance_test(offsets, lbound, window_size, change_point, threshold=1e-15):
    """
    Perform a significance test on a candidate change point using the provided ClaSP.

    Parameters
    ----------
    clasp : ClaSP
        A fitted ClaSP object to use for the significance test.
    change_point : int
        The candidate change point to test.
    threshold : float, optional (default=1e-15)
        The p-value threshold for significance. The default value is 1e-15.

    Returns
    -------
    bool
        True if the change point is significant, False otherwise.

    Notes
    -----
    This method uses a two-sample rank-sum test with a p-value threshold to determine if a candidate change point
    is statistically significant. The test is performed using the classification labels for the candidate change point.
    If the resulting p-value is less than or equal to the specified threshold, the change point is considered significant
    and the method returns True. Otherwise, the change point is considered not significant and the method returns False.

    """
    _, y_pred = cross_val_labels(offsets, change_point-lbound, window_size)
    _, p = ranksums(y_pred[:change_point], y_pred[change_point:])
    return p <= threshold


def score_threshold(profile, change_point, threshold=0.75):
    """
    Returns whether the ClaSP score at the given change point exceeds the specified threshold.

    Parameters
    ----------
    clasp : ClaSP object
        An instance of the ClaSP class that has already been fit to a time series.
    change_point : int
        The index of the change point to be tested.
    threshold : float, optional
        The threshold value to test against the ClaSP score. The default is 0.75.

    Returns
    -------
    bool
        True if the ClaSP score at the given change point is greater than or equal to the specified
        threshold, False otherwise.
    """
    return profile[change_point] >= threshold

def cross_val_labels(offsets, split_idx, window_size):
    """
    Generate predicted and true labels for cross-validation based on nearest neighbour distances.

    Parameters
    ----------
    offsets : ndarray of shape (n_timepoints, k_neighbours)
        The indices of the nearest neighbours for each timepoint in the time series. These indices
        are relative to the start of the time series and should be positive integers.
    split_idx : int
        The index at which to split the time series into two potential segments. This index should be
        less than n_timepoints and greater than window_size.
    window_size : int
        The size of the window used to calculate nearest neighbours.

    Returns
    -------
    y_true : ndarray of shape (n_timepoints,)
        The true labels for each timepoint in the time series.
    y_pred : ndarray of shape (n_timepoints,)
        The predicted labels for each timepoint in the time series.
    """
    n_timepoints, k_neighbours = offsets.shape
    print(n_timepoints)
    print(k_neighbours)
    print(split_idx)

    y_true = np.concatenate((
        np.zeros(split_idx, dtype=np.int64),
        np.ones(n_timepoints - split_idx, dtype=np.int64),
    ))
    print(y_true.shape)

    knn_labels = np.zeros(shape=(k_neighbours, n_timepoints), dtype=np.int64)
    print(knn_labels.shape)

    for i_neighbor in range(k_neighbours):
        neighbours = offsets[:, i_neighbor]
        print(neighbours)
        print(type(neighbours))
        knn_labels[i_neighbor] = y_true[neighbours]

    ones = np.sum(knn_labels, axis=0)
    zeros = k_neighbours - ones
    y_pred = np.asarray(ones > zeros, dtype=np.int64)

    exclusion_zone = np.arange(split_idx - window_size, split_idx)
    y_pred[exclusion_zone] = 1

    return y_true, y_pred