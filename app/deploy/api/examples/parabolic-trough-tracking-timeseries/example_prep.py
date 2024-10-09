#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 14:21:14 2023

@author: bstanisl

Create demo for NREL/SolTrace pull request
"""
import pandas as pd
from st_processing_functions import *

def find_year_month_day(times):
    if times.year[0]==times.year[-1]:
        year = str(times.year[0])
    else:
        year = '*'
    
    if times.month[0]==times.month[-1]:
        month = times.strftime('%m')[0]
    else:
        month = '*'
        
    if times.day[0]==times.day[-1]:
        day = times.strftime('%d')[0]
    else:
        day = '*'
    
    return year,month,day

# prepare demo field data
# copied from run_pysoltrace_iterate.py
tstart = '2023-01-15 16:00:00' # fulldata.index[0] # '2023-02-11 17:00:00'
tend = '2023-01-15 21:00:00' 
times = pd.date_range(tstart, tend, freq='0.5H') # in UTC
field_data_path = '/Volumes/Processed_data/'

year, month, day = find_year_month_day(times)
fileres = '1min' # '1min' or '20Hz'
outres = '0.5H'

field_data = load_field_data(field_data_path, year, month, day, fileres, outres)

# calculate adjusted tilt angle

# select only necessary columns
# selcols = list(field_data.filter(regex='Tilt_adj|trough_angle_').columns)
# selcols = list(field_data.filter(regex='Tilt_adj').columns)
selcols = list(field_data.filter(regex='Tilt$').columns)
selcols.extend(['wspd_3m','wdir_3m',])
field_data = field_data[selcols]

# sample field data at specified times
field_data = field_data.loc[times]

demo_field_data = field_data[['R4_SO_Tilt','R4_DO_Tilt']].copy()
demo_field_data.to_pickle('/Users/bstanisl/OneDrive - NREL/Documents/seto-csp-project/SolTrace/s_SolTrace_gitclone_10_31_23/SolTrace/app/deploy/api/demo_field_data.p')
