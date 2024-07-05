import numpy as np
import pandas as pd
import os

list_mode = []
list_thr = []
list_wd = []
list_res = []

input_dir = '/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_nanopore/claspy/median'
results = os.listdir(input_dir)
        
for f in results:
    filename, extension = os.path.splitext(f)
    if extension == '.npy' and filename.split('_')[-1] == 'maf':
        print(filename)
        file_res = os.path.join(input_dir, f)

        mode = filename.split('_')[0]
        wid_size = filename.split('_')[1]
        thr = filename.split('_')[2]
        
        pred_cp = list(np.unique(np.sort(np.load(file_res))))

        list_mode.append(mode)
        list_thr.append(thr)
        list_wd.append(wid_size)
        list_res.append(pred_cp)
        
df_res = pd.DataFrame({'mode':list_mode,
            'thr':list_thr,
            'ws':list_wd,
            'res':list_res})
print(df_res)

df_res.to_csv(input_dir + '/results_maf_median.tsv', 
            sep= "\t", 
            index = False) 
