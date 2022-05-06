#!/usr/bin/env python3

import numpy as np
import pandas as pd
import re

def import_data(path="../assets/GDSC2_fitted_dose_response_25Feb20.csv"):
    '''Import GSDC data
    
    Params:
    - path (str): path to GDSC csv data

    Output:
    - df_wide (pandas.DataFrame): cell x drug dataframe of lnIC50 values
    - mx_gdsc (numpy.ndarray): 
    - list_drug (list): list of strings of drug names
    - list_cell (list): list of strings of cell names
    '''
    import pandas as pd

    df = pd.read_csv(path)

    ## drug unique contains the drug name and ID concatenated
    df['DRUG_UNIQUE'] = df['DRUG_NAME'] + ' (' + df['DRUG_ID'].apply(str) + ')'

    ## create wide dataframe
    df_wide = df.pivot('CELL_LINE_NAME', 'DRUG_UNIQUE', 'LN_IC50')

    ## sort columns by median value
    sorted_index = df_wide.median().sort_values().index
    df_wide = df_wide[sorted_index]

    ## create matrix of LN_IC50 and list of drugs/cells
    mx_gdsc = df_wide.to_numpy()
    list_drug = list(df_wide.columns)
    list_cell = list(df_wide.index)
    
    return df_wide, mx_gdsc, list_drug, list_cell

def process_data(mx_input, percent_val=0.2, transpose=False, seed=123):
    '''Preprocess data and return indices, masks, and centered data
    
    Params:
    - mx_input (numpy.ndarray): input matrix
    - precent_val (float): percent of training data allocated to validation
    - transpose (bool): if True drug x cell else cell x drug
    - seed (int): set seed for reproducibility
    
    Output:
    - mx_center (numpy.ndarray): matrix of lnIC50 centered by drug (drug x cell if transpose=True, else cell x drug) 
    - mx_train (numpy.ndarray): same as above but using mean from training data and imputing for val/test splits 
    - idx_val (numpy.ndarray): indices of validation data
    - idx_train (numpy.ndarray): indices of training data
    - mask_val (numpy.ndarray): boolean mask of validation data
    - mask_train (numpy.ndarray): boolean mask of training data
    '''
    np.random.seed(seed)
    
    n_cell, n_drug = mx_input.shape
    n_total = mx_input.size
    
    if transpose:
        ## drug x cell
        mx_input = mx_input.T
        
    ## masks for missing, flatten indices, and shuffle
    mask_missing = np.isnan(mx_input)
    idx_nonmissing = np.flatnonzero(~mask_missing)
    idx_shuffle = idx_nonmissing[np.random.permutation(np.arange(len(idx_nonmissing)))]

    ## train/val split
    n_val = int(np.rint(np.sum(~mask_missing) * percent_val))
    n_train = int(np.sum(~mask_missing) - n_val)

    ## find index for train/val
    idx_val = idx_shuffle[:n_val]
    idx_train = idx_shuffle[n_val:(n_val + n_train)]

    ## create mask for train/val
    mask_val = np.array([True if x in idx_val else False for x in range(n_total)])
    mask_train = np.array([True if x in idx_train else False for x in range(n_total)])
    if transpose:
        ## drug x cell
        mask_val = mask_val.reshape((n_drug, n_cell))
        mask_train = mask_train.reshape((n_drug, n_cell)) 
    else:
        ## cell x drug
        mask_val = mask_val.reshape((n_cell, n_drug))
        mask_train = mask_train.reshape((n_cell, n_drug))
        
    ## mean center data by drug for total and training
    ## drug x cell
    if transpose:
        ## for total data
        mean_total = np.nanmean(mx_input, axis=1)[:, np.newaxis]
        mx_center = np.where(mask_missing, mean_total, mx_input)
        mx_center = mx_center - mean_total
        ## for training data
        mean_train = np.ma.mean(np.ma.array(mx_input, mask=~mask_train), axis=1)[:, np.newaxis]
        mx_train = np.where(~mask_train, mean_train, mx_input)
        mx_train = mx_train - np.mean(mx_train, axis = 1)[:, np.newaxis]
    ## cell x drug
    else:
        ## for total data
        mean_total = np.nanmean(mx_input, axis=0)
        mx_center = np.where(mask_missing, mean_total, mx_input)
        mx_center = mx_center - mean_total
        ## for training data
        mean_train = np.ma.mean(np.ma.array(mx_input, mask=~mask_train), axis=0)
        mx_train = np.where(~mask_train, mean_train, mx_input)
        mx_train = mx_train - np.mean(mx_train, axis = 0)
        
    return mx_center, mx_train, idx_val, idx_train, mask_val, mask_train

def process_test(mx_input, percent_test=0.2, percent_val=0.2, transpose=False, seed=123):
    '''Preprocess data and return indices, masks, and centered data
    
    Params:
    - mx_input (numpy.ndarray): input matrix
    - precent_test (float): percent of total non-missing data allocated to testing 
    - precent_val (float): percent of training data allocated to validation
    - transpose (bool): if True drug x cell else cell x drug
    
    Output:
    - mx_center (numpy.ndarray): matrix of lnIC50 centered by drug (drug x cell if transpose=True, else cell x drug) 
    - mx_train (numpy.ndarray): same as above but using mean from training data and imputing for val/test splits 
    - idx_test (numpy.ndarray): indices of test data
    - idx_val (numpy.ndarray): indices of validation data
    - idx_train (numpy.ndarray): indices of training data
    - mask_test (numpy.ndarray): boolean mask of test data
    - mask_val (numpy.ndarray): boolean mask of validation data
    - mask_train (numpy.ndarray): boolean mask of training data
    '''
    np.random.seed(seed)
    
    n_cell, n_drug = mx_input.shape
    n_total = mx_input.size
    
    if transpose:
        ## drug x cell
        mx_input = mx_input.T
        
    ## masks for missing, flatten indices, and shuffle
    mask_missing = np.isnan(mx_input)
    idx_nonmissing = np.flatnonzero(~mask_missing)
    idx_shuffle = idx_nonmissing[np.random.permutation(np.arange(len(idx_nonmissing)))]

    ## train/val/test split
    n_test = int(np.rint(np.sum(~mask_missing) * percent_test))
    n_val = int(np.rint(np.sum(~mask_missing) * (1 - percent_test) * percent_val))
    n_train = int(np.sum(~mask_missing) - n_test - n_val)

    ## find index for train/val/test
    idx_test = idx_shuffle[:n_test]
    idx_val = idx_shuffle[n_test:(n_test + n_val)]
    idx_train = idx_shuffle[(n_test + n_val):(n_test + n_val + n_train)]

    ## create mask for train/val/test
    mask_test = np.array([True if x in idx_test else False for x in range(n_total)])
    mask_val = np.array([True if x in idx_val else False for x in range(n_total)])
    mask_train = np.array([True if x in idx_train else False for x in range(n_total)])
    ## drug x cell
    if transpose:
        mask_test = mask_test.reshape((n_drug, n_cell))
        mask_val = mask_val.reshape((n_drug, n_cell))
        mask_train = mask_train.reshape((n_drug, n_cell)) 
    ## cell x drug
    else:
        mask_test = mask_test.reshape((n_cell, n_drug))
        mask_val = mask_val.reshape((n_cell, n_drug))
        mask_train = mask_train.reshape((n_cell, n_drug))
        
    ## mean center data by drug for total and training
    ## drug x cell
    if transpose:
        ## for total data
        mean_total = np.nanmean(mx_input, axis=1)[:, np.newaxis]
        mx_center = np.where(mask_missing, mean_total, mx_input)
        mx_center = mx_center - mean_total
        ## for training data
        mean_train = np.ma.mean(np.ma.array(mx_input, mask=~mask_train), axis=1)[:, np.newaxis]
        mx_train = np.where(~mask_train, mean_train, mx_input)
        mx_train = mx_train - np.mean(mx_train, axis = 1)[:, np.newaxis]
    ## cell x drug
    else:
        ## for total data
        mean_total = np.nanmean(mx_input, axis=0)
        mx_center = np.where(mask_missing, mean_total, mx_input)
        mx_center = mx_center - mean_total
        ## for training data
        mean_train = np.ma.mean(np.ma.array(mx_input, mask=~mask_train), axis=0)
        mx_train = np.where(~mask_train, mean_train, mx_input)
        mx_train = mx_train - np.mean(mx_train, axis = 0)
        
    return mx_center, mx_train, idx_test, idx_val, idx_train, mask_test, mask_val, mask_train

def load_cell(path='../assets/cell_list.csv'):
    '''Import cell info data
    
    Params:
    - path (str): path to cell info csv data

    Output:
    - df (pandas.DataFrame): pre-processed dataframe containing additional cell line info
    '''

    df = pd.read_csv(path)

    ## only keep GDSC2 data
    df = df[df['Dataset'] == 'GDSC2']

    ## set NCCIT to unclassified instead of NaN
    df['TCGA Classification'][df['Name'] == 'NCCIT'] = 'UNCLASSIFIED'

    return df

def load_drug(drug_list, path='../assets/drug_data.csv'):
    '''Import drug info data
    
    Params:
    - path (str): path to drug info csv data

    Output:
    - df (pandas.DataFrame): dataframe containing additional drug info
    '''

    df_drug_info = pd.read_csv('../assets/drug_data.csv')

    dict_short = {'pathway_name': list(df_drug_info['pathway_name'].unique()),
                  'pathway_short': ['Pro stab/deg', 'PI3K/MTOR', 'DNA rep',
                                    'EGFR', 'Mitosis', 'Cell cycle', 'Other',
                                    'RTK', 'Apoptosis', 'Unclassified', 'MAPK',
                                    'Misc kinase', 'WNT', 'Genome', 'Hist meth',
                                    'Chromatin other', 'Metabolism', 'p53',
                                    'Cytoskeleton', 'Hormone', 'Hist acet', 'IGF1R',
                                    'ABL', 'JNK/p38']}

    df_path = pd.DataFrame(data=dict_short)

    drug_list_sub = [re.sub("\)", "", string) for string in [re.sub(" \(", "\n", string) for string in drug_list]]

    list_replace_to = ['KRAS (G12C) Inhibitor-12\n1855', 'Nutlin-3a (-)\n1047']
    list_replace_from = ['KRAS\nG12C Inhibitor-12\n1855', 'Nutlin-3a\n-\n1047']

    for idx, drug in enumerate(drug_list_sub):
        if drug in list_replace_from:
            drug_list_sub[idx] = list_replace_to[int(np.where(np.array(list_replace_from) == drug)[0])]
            
    drug_list_sub = [string.split('\n')[0] for string in drug_list_sub]

    df = pd.DataFrame({'drug': drug_list_sub,
                       'drug_unique': drug_list})

    df = df.merge(df_drug_info[['drug_name', 'pathway_name']], how='left', left_on='drug', 
                  right_on='drug_name').drop_duplicates().reset_index(drop=True)

    df = df.merge(df_path, how='left', on='pathway_name')

    return df