"""
Basic batch script for pulling data from medical claims files.
"""

import multiprocessing
import argparse
import csv
import math
import os
import pickle
import shutil
from sys import exit

import numpy as np
import pandas as pd


# TODO move this into a config file somewhere

proj_dir = '/project/liu_optum_causal_inference/bp-threshold'

"""diagnosis codes prefixes"""
# heart attack
icd9_mi = ['410']
icd10_mi = ['I21', 'I22']

# hypertension
icd9_ht = ['401']
icd10_ht = ['I10']

# diabetes
icd9_dm = ['250']
icd10_dm = ['E10', 'E11', 'E12', 'E13', 'E14']

# stroke
icd9_st = ['433', '434', '435']
icd10_st = ['I63', 'I64', 'G45']

# retinopathy
icd9_rt = ['3620']
icd10_rt = ['E103', 'E113']


"""Lab reading prefixes"""
# BP readings
loinc_bp = ['55284-4', '8480-6', '8462-4', '35094-2', '85354-9', '8478-0']
cpt_bp = ['3074F', '3075F', '3077F', '3078F', '3079F', '3080F']

# A1C readings
loinc_a1c = ['4548-4', '4549-2', '17856-6', '59261-8', '62388-4', '41995-2']

#with open("{}/pkl_data/diag_cols.list".format(proj_dir), "rb") as diag_f:
#    diag_cols = pickle.load(diag_f)

#with open("pkl_data/proc_cols.list", "rb") as proc_f:
#    proc_cols = pickle.load(proc_f)

proc_cols = ['proc_cd']
keep_cols = ['patid', 'fst_dt']


def diag_extract(in_file, out_path, diag_dict, chunksize, test=False, bp=False):
    """Extracts diagnosis information from medical claim files (m).

    Filters by the ICD9/10 prefixes provided in diag_dict.

    Parameters
    ----------
    in_file : str
        input file path
    out_path : str
        output file path
    diag_dict: dict
        dictionary with (out_name, [icd prefixes]) k,v pairs. eg "ht": ["401", "I10"]
    chunksize : int
        number of entries in the raw file to read per iteration, don't set too large
    test : bool
        enables debugging print statements
    """

    columns = keep_cols + diag_cols

    diag_dfs = {diag: pd.DataFrame() for diag in diag_dict.keys()}

    for chunk in pd.read_stata(in_file, chunksize=chunksize, columns=columns):
        if test: print(chunk.shape[0])
        diag_masks = {diag: [False] * chunk.shape[0] for diag in diag_dict.keys()}

        # check every diagnosis column whether or not the target diagnosis diag is present
        for diag, icd_codes in diag_dict.items():
            for col in diag_cols:
                # bitwise OR every col mask with the diagnosis
                for code in icd_codes:
                    diag_masks[diag] = diag_masks[diag] | chunk[col].str.startswith(code)
            if diag_masks[diag].sum() > 0:
                diag_dfs[diag] = diag_dfs[diag].append(chunk[diag_masks[diag]])
                
        if test:
            for diag, df in diag_dfs.items():
                print(diag)
                print(df.shape)
                print(df.head())
            break

    for diag, df in diag_dfs.items():
        pickle.dump(df, open("{}_{}.df".format(out_path, diag), "wb"), -1)

        
def m_extract(in_file, out_name, chunksize, test=False, bp=False):
    """Old m table processing function."""
    #columns = keep_cols + diag_cols + proc_cols
    columns = keep_cols + diag_cols

    if bp:
        columns = keep_cols + proc_cols
    
    
    mi_df = pd.DataFrame()
    ht_df = pd.DataFrame()
    bp_df = pd.DataFrame()
    

    for chunk in pd.read_stata(in_file, chunksize=chunksize, columns=columns):
        print(chunk.shape[0])
        mi_mask = [False] * chunk.shape[0]
        ht_mask = [False] * chunk.shape[0]
        for col in diag_cols:
            mi_mask = mi_mask | chunk[col].str.startswith(icd9_mi) | chunk[col].str.startswith(icd10_mi[0]) | chunk[col].str.startswith(icd10_mi[1])
            ht_mask = ht_mask | chunk[col].str.startswith(icd9_ht) | chunk[col].str.startswith(icd10_ht) 

        if mi_mask.sum() > 0:
            mi_df = mi_df.append(chunk[mi_mask][keep_cols])

        if ht_mask.sum() > 0:
            ht_df = ht_df.append(chunk[ht_mask][keep_cols])

        # bp_chunk = chunk.loc[chunk['proc_cd'].isin(cpt_bp)]
        # if bp_chunk['proc_cd'].count() > 0:
        #     bp_df = bp_df.append(bp_chunk[keep_cols + ['proc_cd']])
                
        if test:
            print(mi_df.shape)
            print(mi_df.head())

            print(ht_df.shape)
            print(ht_df.head())
            break
            # print(bp_df.shape)
            # print(bp_df.head())

    pickle.dump(mi_df, open(out_name + "_mi.df", "wb"), -1)
    pickle.dump(ht_df, open(out_name + "_ht.df", "wb"), -1)
    # pickle.dump(bp_df, open(out_name + "_bp.df", "wb"), -1)


# TODO parameterize
def lr_extract(in_file, out_name, chunksize, test=False):
    """Extracts lab report information from lab report files (lr).

    Currently extracts BP lab reports and A1C lab reports.
    """
    lab_cols = ['patid', 'pat_planid', 'fst_dt', 'loinc_cd', 'rslt_nbr', 'rslt_txt', 'rslt_unit_nm']


    bp_df = pd.DataFrame()
    a1c_df = pd.DataFrame()


    for chunk in pd.read_stata(in_file, chunksize=chunksize, columns=lab_cols):
        print(chunk.shape[0])
        
        bp_chunk = chunk.loc[chunk['loinc_cd'].isin(loinc_bp)]
        if bp_chunk['loinc_cd'].count() > 0:
            bp_df = bp_df.append(bp_chunk)

        a1c_chunk = chunk.loc[chunk['loinc_cd'].isin(loinc_a1c)]
        if a1c_chunk['loinc_cd'].count() > 0:
            a1c_df = a1c_df.append(a1c_chunk)

            
        if test:
            print(bp_df.shape)
            print(bp_df.head())

            print(a1c_df.shape)
            print(a1c_df.head())
            break
            
    pickle.dump(bp_df, open(out_name + "_lr_bp.df", "wb"), -1)
    pickle.dump(a1c_df, open(out_name + "_lr_a1c.df", "wb"), -1)


def rx_extract(in_file, out_path, out_name, ndcs, chunksize, test=False):
    """Extracts rx pharmacy information from (r) files.

    Currently filters by the provided list of ndc codes.
    
    Parameters
    ----------
    in_file : str
        input file path
    out_path : str
        output file path
    out_name : str
        output file name, usually something short (eg dm for diabetes mellitus)
    ndcs : list
        list of ndcs to filter by
    chunksize : int
        number of entries in the raw file to read per iteration, don't set too large
    test : bool
        enables debugging print statements

    Returns
    -------
    None 
        output is dumped to pickle files
    """
    
    # columns to extract
    rx_cols = ['patid', 'ndc', 'npi', 'chk_dt', 'fill_dt', 'fst_fill']
    ndc_col = 'ndc'
    
    # dumb discrepancy where the 2001 q1 columns are capitalized
    if '2001q1' in out_path:
        rx_cols = ['Patid', 'Ndc', 'Npi', 'Chk_Dt', 'Fill_Dt', 'Fst_Fill']
        ndc_col = 'Ndc'
        
    rx_df = pd.DataFrame()

    for chunk in pd.read_stata(in_file, chunksize=chunksize, columns=rx_cols):
        if test: print(chunk.shape[0])
        
        rx_chunk = chunk.loc[chunk[ndc_col].isin(ndcs)]
        rx_df = rx_df.append(rx_chunk)
            
        if test:
            print(rx_df.shape)
            print(rx_df.head())
            break
            
    pickle.dump(rx_df, open("{}_rx_{}.df".format(out_path, out_name), "wb"), -1)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract data from Optum raw files and dump to DataFrames")
    parser.add_argument('data_dir', help='directory with all Optum data')
    parser.add_argument('yr', help='the year to target')
    parser.add_argument('q', help='the quarter to target')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('table_type', choices=['m', 'lr', 'r'], help='Optum table type to target: medical (m), lab reports (lr), prescriptions (r)') 
    parser.add_argument('chunksize', type=int, help='number of rows to read per chunk')
    parser.add_argument('--test', action='store_true', help='whether to make a test run of the data extraction')
    parser.add_argument('--m_dm_outcome', action='store_true', help='perform diabetes med outcome extraction')
    parser.add_argument('--rx_dm', action='store_true', help='perform diabetes rx extraction')

    
    args = parser.parse_args()
    
    name = "ses_{}{}q{}".format(args.table_type, args.yr, args.q)
    i_file = args.data_dir + name + ".dta"
    o_file = args.out_dir + name 
            
    # for q in range(1, 5):
    #     if args.table_type == 'lr':
    #         names.append("ses_lr{}q{}".format(args.yr, q))
    #     elif args.table_type == 'r':
    #         names.append("ses_r{}q{}".format(args.yr, q))
    #     elif args.table_type == 'm': 
    #         names.append("ses_m{}q{}".format(args.yr, q))
    # i_files = [args.data_dir + name + ".dta" for name in names]
    # o_files = [args.out_dir + name for name in names]
    # extract_args = [(i_file, o_file, args.chunksize, args.test, args.bp) for i_file, o_file in zip(i_files, o_files)]

    print(i_file)
    print(o_file)

    if args.table_type == 'm':
        if args.m_dm_outcome:            
            extract_fun = diag_extract
            diag_dict = {
                "stroke": icd9_st + icd10_st,
                "diab": icd9_dm + icd10_dm,
                "retin": icd9_rt + icd10_rt,
                "h_atk": icd9_mi + icd10_mi
            }
            if args.test: print(diag_dict)

            extract_fun(i_file, o_file, diag_dict, args.chunksize, args.test)
            
    elif args.table_type == 'lr':
        extract_fun = lr_extract
        # TODO fix labs flag

    elif args.table_type == 'r':
        extract_fun = rx_extract

        if args.rx_dm:
            ndcs = pickle.load(open("{}/pkl_data/dm_ndc.list".format(proj_dir), "rb"))
            o_name = 'dm'
            
        extract_args = [(i_file, o_file, o_name, ndcs, args.chunksize, args.test) for i_file, o_file in zip(i_files, o_files)]

        for i_f, o_f, o_n, ndcs, cs, test in extract_args:
            rx_extract(i_f, o_f, o_n, ndcs, cs, test)
    
