import pandas as pd
import numpy as np

# specific function to call if input file is a .tsv
def get_data_tsv(data: str):
    # return dict to keep track of variable names
    ret_dict = {}
    snv = pd.read_csv(data, sep = '\t')
    snv.fillna(value=0, axis = 0, inplace = True)
    
    ret_dict["vaf"] = np.array(snv.gt_AF)
    ret_dict["baf"] = np.array(snv.median_baf)
    ret_dict["dr"] = np.array(snv.median_dr)
    ret_dict["maf"] = np.array(snv.median_meth)
    
    return ret_dict

# specific funtion to call if a function is a .csv
def get_data_csv(data: str):
    # return dict to keep track of variable names
    ret_dict = {}

    snv = pd.read_csv(data)
    snv.fillna(value='normal', axis = 0, inplace = True)
    snv.sort_values(by = ['pos'], inplace = True)
    snv['id'] = [i for i in range(snv.shape[0])]
                    
    bps_max = list(snv.groupby(['cna_id']).max(['pos']).id)
    bps_min = list(snv.groupby(['cna_id']).min(['pos']).id)               
    vaf_bps = snv.groupby('cna_id')
                    
    ret_dict["vaf"] = np.array(snv.vaf)
    ret_dict["baf"] = np.array(snv.median_baf)
    ret_dict["maf"] = np.array(snv.median_meth)
    ret_dict["dr"] = np.array(snv.median_dr)
    ret_dict["bps"] = np.array(bps_max + bps_min) 
    
    return ret_dict
