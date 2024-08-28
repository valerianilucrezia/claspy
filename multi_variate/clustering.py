from sklearn.metrics import silhouette_samples
import numpy as np


def cluster_segments(ts: np.ndarray, cp: np.ndarray) -> np.ndarray:
    """Find binary segmented time series clusters of mutated vs non-mutated regions.

    Args:
        ts (np.ndarray): Time series of shape (t,2).
        cp (np.ndarray): List of change point locations for ts.

    Returns
        (np.ndarray): List of size ts of labeled time points. 0 for cluster A, 1 for cluster B.
    """

    cluster_a_cp = sorted(cp)
    cluster_b_cp = sorted(cp)
    cluster_a_cp.insert(0, 0)

    # good because cp lists will never both be of even length
    if len(cluster_a_cp) % 2 == 1:
        cluster_a_cp.append(-1)
    elif len(cluster_b_cp) % 2 ==1:
        cluster_b_cp.append(-1)

    # TODO what cluster is a changepoint in?
    cluster_b_idx = list(zip(*[iter(cluster_b_cp)]*2))

    labels = np.zeros(shape=ts.shape[0])

    for idx in cluster_b_idx:
        if idx[1] != -1:
            labels[idx[0]:idx[1]] = 1
        else:
            labels[idx[0]:] = 1

    return labels


def calculate_cluster_scores(ts:list | np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Function to return the Sihlouette Score for each time point, which denotes goodness of cluster fit.

    Args:
        ts (list | np.ndarray): Time series of shape (t,2).
        labels (np.ndarray): Cluster assignments for each time point.

    Returns:
        (np.ndarray): Sihlouette Score for each time point.
    """
    
    scores = silhouette_samples(X=np.transpose(ts), labels=labels)
    return scores


def construct_confusion_matrix(labels: np.ndarray, scores: np.ndarray) -> np.ndarray:
    """Function to construct a 'confusion matrix' from Sihlouette Scores

    Args:
        labels (np.ndarray): Labeled time points. 0 for cluster A, 1 for cluster B
        
        scores (np.ndarray): Score between -1 and 1 for cluster assignment fit. > 0 is good, < 0 is bad.

    Returns:
        (np.ndarray) 2x2 array with [0,0]=correct A classifications, [0,1]=incorrect A classifications, [1,0]=incorrect B classifications, [1,1]=correct B classifications.
    """
    aa = 0 
    ab = 0
    ba = 0
    bb = 0
    for i in range(len(scores)-1):
        #incorrect assignment
        if scores[i] <= 0:
            if labels[i] == 0:
                # assigned a but should be b
                ab += 1
            else:
                # assigned b but should be a
                ba += 1
        #correct assignment
        else:
            if labels[i] == 0:
                aa += 1
            else:
                bb += 1
    
    # return confusion matrix
    return np.array(((aa, ab), (ba, bb)))
