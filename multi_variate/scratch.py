from tests.sim_data import sim_data
import numpy as np
size=100
profile, change_points = sim_data(base_freq_range=(0.25, 0.55), mutant_freq_range=(0.9, 0.95), size_segment=10, size_ts=size, prob=(0.7, 0.3))
print(profile)
print(change_points)


np.zeros(shape=(40,))
np.ones(shape=(20,))
np.zeros(shape=(10,))
np.ones(shape=(10,))
np.zeros(shape=(10,))
np.ones(shape=(10,))

labels = np.concatenate((np.zeros(shape=(40,)),
                         np.ones(shape=(20,)),
                         np.zeros(shape=(10,)),
                         np.ones(shape=(10,)),
                         np.zeros(shape=(20,)),
                         np.ones(shape=(10,))), axis=None)