import numpy as np
import pandas as pd
import os

base = '/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_races/claspy'

list_sim = []
list_comb = []
list_mode = []
list_thr = []
list_wd = []
list_res = []

for s in range(3,14):
    sim = f'sim_{s}'
    path = os.path.join(base, sim)
    combinations = os.listdir(path)
    
    for comb in combinations:
        if comb[0] == 'c':
            path_comb = os.path.join(path, comb)
            results = os.listdir(path_comb) 
            for f in results:
                filename, extension = os.path.splitext(f)
                if extension == '.npy':
                    file_res = os.path.join(path_comb, f)
        
                    mode = filename.split('_')[0]
                    wid_size = filename.split('_')[1]
                    thr = filename.split('_')[2]
                    
                    pred_cp = list(np.unique(np.sort(np.load(file_res))))

                    list_sim.append(sim)
                    list_comb.append(comb)
                    list_mode.append(mode)
                    list_thr.append(thr)
                    list_wd.append(wid_size)
                    list_res.append(pred_cp)
                    
df_res = pd.DataFrame({'sim_id':list_sim,
                        'comb':list_comb,
                        'mode':list_mode,
                        'thr':list_thr,
                        'ws':list_wd,
                        'res':list_res})
print(df_res)

df_res.to_csv('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_races/claspy/results.tsv', 
              sep= "\t", 
              index = False) 


list_sim = []
list_bp = []
for s in range(3,14):
    sim = f'sim_{s}'
    path = os.path.join(base, sim)
    real_bp = list(np.load(path + '/bp.npy'))
    
    list_sim.append(sim)
    list_bp.append(real_bp)
    
df_real = pd.DataFrame({'sim_id':list_sim,
                        'real_bp':list_bp})
df_real.to_csv('/orfeo/LTS/LADE/LT_storage/lvaleriani/CNA/segmentation/res_races/claspy/bps.tsv', 
              sep= "\t", 
              index = False)  
                    
                    
                    