"""
This script is part of the TIMESAT toolkit for the RAMONA project.

Author: Zhanzhang Cai

This Python script is designed to perform a variety of tasks related to time-series analysis of satellite sensor data, in the context of the RAMONA project.

Instructions for running the code:
1. Ensure that you have the necessary packages installed and the data available in the correct directories.
2. In the terminal, navigate to the directory containing this script.
3. Enter 'python run.py -help'.
4. Follow any prompts that appear in the terminal.

Please ensure you have read the associated documentation and understand the functions within this code before running it.

Please modify L20 in settings_test_evi2.json for the path of output.
For test:
Enter 'python run.py -i settings_test_evi2.json' in your terminal

For any further questions or assistance, contact Zhanzhang Cai.
"""
# import sys
# sys.path.append("/Users/zhanzhangcai/Documents/GitHub/TIMESAT_python")  # Replace with the actual path
# from memory_profiler import profile
import json
import math
import os
import numpy as np
import timesat
import time
import copy
import re
from datetime import datetime, timedelta, date
import pandas as pd
import matplotlib.pyplot as plt
import csv, gc
import glob

def assign_qa_weight(p_a, qa):
    p_a = np.asarray(p_a)
    
    # If p_a is empty, return the default qa_out.
    if p_a.size == 0:
        return qa

    # Create a new output array with default value 0.
    qa_out = np.zeros_like(qa, dtype=float)
    
    if p_a.shape[1] == 2:
        # Mapping: [qa_value, weight]
        for rule in p_a:
            qa_value, weight = rule
            mask = (qa == qa_value)
            qa_out[mask] = weight
    elif p_a.shape[1] == 3:
        # Mapping: [min_value, max_value, weight]
        for rule in p_a:
            min_val, max_val, weight = rule
            mask = (qa >= min_val) & (qa <= max_val)
            qa_out[mask] = weight
    else:
        raise ValueError("p_a must have either 2 or 3 columns.")
    
    return qa_out


def yydoy_float_to_datetime(sosd_day):
    """
    Convert a YYDOY float (e.g., 24123.5 -> 2024-05-02 12:00)
    into a Python datetime object.
    Handles fractional days and two-digit years.
    """

    yy = int(sosd_day // 1000)
    doy_float = sosd_day % 1000  # retains decimals (e.g., 123.5)


    # Interpret 2-digit year: 00–69 => 2000–2069, 70–99 => 1970–1999
    year = 2000 + yy

    doy_int = int(doy_float)

    sos_date = datetime(year, 1, 1) + timedelta(days=doy_int - 1)

    return sos_date



def assign_qa_weight(p_a, qa):
    p_a = np.asarray(p_a)
    
    # If p_a is empty, return the default qa_out.
    if p_a.size == 0:
        return qa

    # Create a new output array with default value 0.
    qa_out = np.zeros_like(qa, dtype=float)
    
    if p_a.shape[1] == 2:
        # Mapping: [qa_value, weight]
        for rule in p_a:
            qa_value, weight = rule
            mask = (qa == qa_value)
            qa_out[mask] = weight
    elif p_a.shape[1] == 3:
        # Mapping: [min_value, max_value, weight]
        for rule in p_a:
            min_val, max_val, weight = rule
            mask = (qa >= min_val) & (qa <= max_val)
            qa_out[mask] = weight
    else:
        raise ValueError("p_a must have either 2 or 3 columns.")
    
    return qa_out


def read_table_data(df, p_a):

    # Identify the first column dynamically
    first_col_name = df.columns[0]

    # Convert the first column to datetime
    dates = pd.to_datetime(df[first_col_name])

    # Convert dates to the desired format yyyyDOY
    tv_yyyydoy = dates.dt.strftime('%Y%j').tolist()

    # Convert the first column to the desired format (assuming it contains dates)
    tv_yyyymmdd = dates.dt.strftime('%Y%m%d').tolist()

    # Get the number of rows
    npt = len(df)

    col = 1

    # Initialize 3D arrays
    ym3 = np.zeros((col, 1, npt), dtype=np.float32)
    wm3 = np.ones((col, 1, npt), dtype=np.float32)
    lc2 = np.ones((1, 1), dtype=np.float32)

    # Save data from the second column to the end into ym3
    # Transpose the slice to align with the shape (col, 1, npt)
    ym3[0, 0, :] = df.iloc[:, 1].T.values.astype(np.float32)
    try:
        wm3[0, 0, :] = df.iloc[:, 2].values.astype(np.float32)
    except:
        pass


    # Convert lists to numpy arrays
    tv_yyyymmdd = np.array(tv_yyyymmdd, dtype=int)
    tv_yyyydoy = np.array(tv_yyyydoy, dtype=int)

    # Extract year and calculate the number of unique years from YYYYMMDD format
    yv = tv_yyyymmdd // 10000  # Extract year from YYYYMMDD (integer division by 10000)
    yrstart = np.min(yv)  # First year
    yrend = np.max(yv)  # First year
    nyear = np.max(yv) - yrstart + 1  # Number of unique years (inclusive)

    # Compute qa_weight using the 2-column mapping
    wm3 = assign_qa_weight(p_a, wm3)

    return tv_yyyydoy, tv_yyyymmdd, ym3, wm3, lc2, nyear, yrstart, yrend, npt, col


def yydoy_float_to_datetime(sosd_day):
    """
    Convert a YYDOY float (e.g., 24123.5 -> 2024-05-02 12:00)
    into a Python datetime object.
    Handles fractional days and two-digit years.
    """

    yy = int(sosd_day // 1000)
    doy_float = sosd_day % 1000  # retains decimals (e.g., 123.5)


    # Interpret 2-digit year: 00–69 => 2000–2069, 70–99 => 1970–1999
    year = 2000 + yy

    doy_int = int(doy_float)

    sos_date = datetime(year, 1, 1) + timedelta(days=doy_int - 1)

    return sos_date


def _ts_run_(indata,run_vpp):

    viname = [c for c in indata.columns if c != "t"][0]

    if viname == 'GPP':
        fluxgpp = 1
    else:
        fluxgpp = 0


    p_ignoreday = 366
    # NDVI,NIRv,EVI2,PPI,LAI,FAPAR,GPP_m,GPP
    if fluxgpp or 'GPP_m' in indata:
        p_ylu = [0., 9999.]
    elif 'LAI' in indata:
        p_ylu = [0, 10]
    elif 'NDVI' in indata or 'NIRv' in indata:
        p_ylu = [0, 1]
    elif 'PPI' in indata:
        p_ylu = [0, 5]
    elif 'EVI2' in indata:
        p_ylu = [0, 1]
    elif 'FAPAR' in indata:
        p_ylu = [0, 1]

    p_a = [[-10000., 10000., 1.], [-10000., 10000., 1.], [-10000., 10000., 1.]]
    p_outststep = 1
    p_nodata = -9999
    p_davailwin = 45
    p_outlier = 0


    p_hrvppformat = 1
    p_printflag = 0
    p_nclasses = 1

    # Initialize TIMESAT parameter vectors (length 255)
    landuse = np.zeros(255, dtype='uint8')
    p_fitmethod = np.zeros(255, dtype='uint8')
    p_smooth_vec = np.zeros(255, dtype='double')
    p_nenvi_vec = np.ones(255, dtype='uint8')
    p_wfactnum_vec = np.ones(255, dtype='double')
    p_startmethod_vec = np.ones(255, dtype='uint8')
    p_startcutoff_vec = np.full((255, 2), 0.5, order='F', dtype='double')
    p_low_percentile_vec = np.full(255, 0.05, dtype='double')
    p_fillbase_vec = np.ones(255, dtype='uint8')
    p_seasonmethod_vec = np.zeros(255, dtype='uint8')
    p_seapar_vec = np.zeros(255, dtype='double')

    landuse[0] = 1

    timevector, tv_yyyymmdd, vi, qa, lc, yr, yrstart, yrend, npt, col = read_table_data(indata, p_a)

    p_outindex = np.arange(1, yr * 365 + 1)[::p_outststep]
    p_outindex_num = len(p_outindex)

    daily_timestep = []
    for year in range(yrstart, yrstart + yr):
        for day in range(365):
            daily_timestep.append(datetime(year, 1, 1) + timedelta(days=day))
    daily_timestep = daily_timestep[::p_outststep]

    # --- Prepare static tables (no repeated concat) ---
    seasons = ["s1", "s2"]
    vppname = ["SOSD", "SOSV", "LSLOPE", "EOSD", "EOSV", "RSLOPE", "LENGTH",
               "MINV", "MAXD", "MAXV", "AMPL", "TPROD", "SPROD"]
    vpplist = [
        f"{year}_{season}_{name}"
        for year in range(yrstart, yrend + 1)
        for season in seasons
        for name in vppname
    ]

    # Parameter grids
    if fluxgpp:
        p_fitmethod_options = [2]
        p_smooth_options_method_list = [5000]
        p_seasonmethod_options_list = [1]
        p_seapar_options_list = [0.5]
        p_startcutoff_options = [0.4]
        p_nenvi_vec[0] = 3
        p_wfactnum_vec[0] = 10
    else:
        p_fitmethod_options = [1, 2]
        p_smooth_options_method_list = [100, 300, 1000, 3000, 10000]
        p_seasonmethod_options_list = [1, 2]
        p_seapar_options_list = [0.0, 0.25, 0.5, 0.75, 1.0]
        if run_vpp == 1:
            p_startcutoff_options = np.arange(0.05, 0.55, 0.05)
        elif run_vpp == 0:
            p_startcutoff_options = [0.4]

    # Collect columns in memory
    vpp_columns = []        # list of 1D arrays length len(vpplist)
    yfit_columns = []       # list of 1D arrays length len(daily_timestep)
    settings_rows = []      # list of dicts (one per setting)

    counter_st = 0

    for fitmethod_value in p_fitmethod_options:
        p_seasonmethod_options = p_seasonmethod_options_list
        p_seapar_options = p_seapar_options_list
        if fitmethod_value == 1:
            method_name = 'DL'
            p_smooth_options = [p_smooth_options_method_list[0]]   
        elif fitmethod_value == 2:
            method_name = 'SP'
            p_smooth_options = p_smooth_options_method_list
            if run_vpp == 0:
                p_seasonmethod_options = [p_seasonmethod_options_list[0]]
                p_seapar_options = [p_seapar_options_list[0]]


        p_fitmethod.fill(fitmethod_value)

        for p_smooth_val in p_smooth_options:
            p_smooth_vec.fill(p_smooth_val)

            for seasonmethod_value in p_seasonmethod_options:
                p_seasonmethod_vec.fill(seasonmethod_value)

                for seapar_value in p_seapar_options:
                    p_seapar_vec.fill(seapar_value)

                    # Record settings row
                    counter_st += 1

                    counter_vpp = 0

                    for soscutoff_value in p_startcutoff_options:
                        p_startcutoff_vec[0, 0] = soscutoff_value

                        for eoscutoff_value in p_startcutoff_options:
                            p_startcutoff_vec[1, 0] = eoscutoff_value

                            counter_vpp += 1

                            if run_vpp == 1:
                                settings_id = viname + f"-{counter_st}-{counter_vpp}"
                            else:
                                settings_id = viname + f"-{counter_st}"

                            settings_rows.append({
                                "settings_id": settings_id,
                                "fitmethod": int(fitmethod_value),
                                "method_name": method_name,
                                "smooth": float(p_smooth_val),
                                "seasonmethod": int(seasonmethod_value),
                                "seapar": float(seapar_value),
                                "sos_cutoff": float(soscutoff_value),
                                "eos_cutoff": float(eoscutoff_value)
                            })

                            if p_printflag == 1:
                                print(counter_st)
                                print(counter_vpp)

                            if 0:
                                print(f"method_name: {method_name}")
                                print(f"p_smooth_val: {p_smooth_val}")
                                print(f"seasonmethod_value: {seasonmethod_value}")
                                print(f"seapar_value: {seapar_value}")
                                print(f"soscutoff_value: {soscutoff_value}")
                                print(f"eoscutoff_value: {eoscutoff_value}")
                                print("\n--- DEBUG: shapes before timesat.tsf2py() ---")
                                def safe_shape(x):
                                    try:
                                        return np.shape(x)
                                    except Exception:
                                        return type(x)

                                debug_vars = {
                                    "yr": yr,
                                    "vi": safe_shape(vi),
                                    "qa": safe_shape(qa),
                                    "timevector": safe_shape(timevector),
                                    "lc": safe_shape(lc),
                                    "p_nclasses": p_nclasses,
                                    "landuse": safe_shape(landuse),
                                    "p_outindex": safe_shape(p_outindex),
                                    "p_ignoreday": p_ignoreday,
                                    "p_ylu": safe_shape(p_ylu),
                                    "p_fitmethod": safe_shape(p_fitmethod),
                                    "p_smooth_vec": safe_shape(p_smooth_vec),
                                    "p_nodata": p_nodata,
                                    "p_outlier": p_outlier,
                                    "p_davailwin": p_davailwin,
                                    "p_nenvi_vec": safe_shape(p_nenvi_vec),
                                    "p_wfactnum_vec": safe_shape(p_wfactnum_vec),
                                    "p_startmethod_vec": safe_shape(p_startmethod_vec),
                                    "p_startcutoff_vec": safe_shape(p_startcutoff_vec),
                                    "p_low_percentile_vec": safe_shape(p_low_percentile_vec),
                                    "p_fillbase_vec": safe_shape(p_fillbase_vec),
                                    "p_hrvppformat": p_hrvppformat,
                                    "p_seasonmethod_vec": safe_shape(p_seasonmethod_vec),
                                    "p_seapar_vec": safe_shape(p_seapar_vec),
                                    "col": col,
                                    "npt": npt,
                                    "p_outindex_num": p_outindex_num,
                                }
                                for k, v in debug_vars.items():
                                    print(f"{k:25s} : {v}")
                                print("--- END DEBUG ---\n")
                                exit()

                            if 1:
                                # Run TIMESAT
                                vpp, vppqa, nseason, yfit, yfitqa, seasonfit, tseq = timesat.tsfprocess(
                                    yr, vi, qa, timevector, lc, p_nclasses, landuse, p_outindex,
                                    p_ignoreday, p_ylu, p_printflag, p_fitmethod, p_smooth_vec, p_nodata,
                                    p_davailwin, p_outlier, p_nenvi_vec, p_wfactnum_vec, p_startmethod_vec,
                                    p_startcutoff_vec, p_low_percentile_vec, p_fillbase_vec, p_hrvppformat,
                                    p_seasonmethod_vec, p_seapar_vec, run_vpp
                                )

                                # Reshape to 1D columns
                                vpp = np.squeeze(np.moveaxis(vpp, -1, 0), axis=-1)

                                if fluxgpp:
                                    vppqa = np.squeeze(np.moveaxis(vppqa, -1, 0), axis=-1)
                                    # print(vpp.shape)
                                    # print(vppqa.shape)
                                    # print(vpp)
                                    # print(vppqa)
                                    n = vppqa.shape[0]
                                    for i in range(n):
                                        if vppqa[i] < 3:
                                            vpp[i*13:(i+1)*13] = p_nodata
                                    # print(vpp)
                                    

                                # Save columns (defer DataFrame creation)
                                vpp_columns.append((settings_id, np.asarray(vpp, dtype=np.float32)))
                                
                                if fluxgpp or run_vpp == 0:
                                    yfit = np.squeeze(np.moveaxis(yfit, -1, 0).astype('float32'), axis=-1)
                                    yfit_columns.append((settings_id, np.asarray(yfit, dtype=np.float32)))
                                    if fluxgpp:
                                        yfitqa = np.squeeze(np.moveaxis(yfitqa, -1, 0).astype('uint8'), axis=-1)
                                        yfit_columns.append((settings_id + '_qa', np.asarray(yfitqa, dtype=np.uint8)))


    # --- Build DataFrames once (fast path) ---
    # Settings table
    settings_df = pd.DataFrame(settings_rows)

    yfit_data = 0
    vpp_data = 0

    if run_vpp == 1:
        # settings_df.to_csv('settings_vpp.csv', index=False)
        # VPP table
        vpp_array = np.column_stack([col_vals for _, col_vals in vpp_columns]) if vpp_columns else np.empty((len(vpplist), 0))
        vpp_cols_names = [sid for sid, _ in vpp_columns]
        vpp_data = pd.DataFrame({"id": vpplist})
        if vpp_array.size > 0:
            vpp_df = pd.DataFrame(vpp_array, columns=vpp_cols_names)
            vpp_data = pd.concat([vpp_data, vpp_df], axis=1)
        vpp_data.replace(-9999, np.nan, inplace=True)
        # vpp_data.to_csv(output_vpp_file_name, index=False)
    elif run_vpp == 0 or fluxgpp:
        # settings_df.to_csv('settings_st.csv', index=False)
        # YFIT table (GPP case only)
        yfit_array = np.column_stack([col_vals for _, col_vals in yfit_columns])
        yfit_cols_names = [sid for sid, _ in yfit_columns]
        yfit_data = pd.DataFrame({"t": daily_timestep})
        yfit_df = pd.DataFrame(yfit_array, columns=yfit_cols_names)
        yfit_data = pd.concat([yfit_data, yfit_df], axis=1)
        # yfit_data.to_csv(output_yfit_file_name, index=False)
        
    return yfit_data, vpp_data, settings_df
