from sklearn.metrics import silhouette_samples
import numpy as np

"""Function to find binary segmented time series clusters and return the Sihlouette Score for each time point, which denotes goodness of cluster fit.
Input:
- ts: time series of shape (t,2)
- cp: list of change point locations on t
Output:
- labels: list of size ts of labeled time points. 0 for cluster_a, 1 for cluster_b
- scores: score between -1 and 1 for cluster assignment fit. > 0 is good, < 0 is bad"""
def calculate_cluster_score(ts: list, cp: list):
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

    labels = np.zeros(ts.shape[1])

    for idx in cluster_b_idx:
        if idx[1] != -1:
            labels[idx[0]:idx[1]] = 1
        else:
            labels[idx[0]:] = 1
    print(ts.shape)
    print(labels)
    # FIXME these are transposed for some reason
    # cluster_a = np.array((ts[:,0][labels == 0], ts[:,1][labels == 0]))
    # cluster_b = np.array((ts[:,0][labels == 1], ts[:,1][labels == 1]))

    scores = silhouette_samples(X=np.transpose(ts), labels=labels)

    return labels, scores


"""Function to construct a 'confusion matrix' from Sihlouette Scores
Input:
- labels: list of size ts of labeled time points. 0 for cluster_a, 1 for cluster_b
- scores: score between -1 and 1 for cluster assignment fit. > 0 is good, < 0 is bad
Output:
- confusion_matrix: matrix of size (2,2) that says which timepoints were correctly assigned to the correct segment or not:
    Format: (
            (Truly A, Assigned A but truly B),
            (Assigned B but truly A, truly B)
            )
"""
def construct_confusion_matrix(labels, scores):
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

# test
baseline = np.random.uniform(0.13, 0.17, 20)
mutant = np.random.uniform(0.53, 0.57, 20)
profile = np.concatenate((baseline, mutant, baseline, mutant, baseline), axis=None)
ts = np.array((range(100),profile))
cp = list(range(21, 100, 20))
print(cp)
labels, scores = calculate_cluster_score(ts, cp)
print(cp)
print(scores)

matrix = construct_confusion_matrix(labels, scores)
print(matrix)

