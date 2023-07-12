#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 08:23:11 2023

@author: bstanisl
"""

# convert pickle files to .csv
import pickle
import pandas as pd

# stats_pickle = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats.p','rb'))
# median_pickle = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/median_diurnal_cycle_allrows.p','rb'))

stats_pickle = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats_no_adj.p','rb'))
median_pickle = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/median_diurnal_cycle_allrows_no_adj.p','rb'))
intercept_pickle = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/stats_1E+05hits_realistic_optics.p', 'rb'))[1]
#%%
stats_pickle2 = stats_pickle.loc[:, ['absmean', 'absstd', 'absmax', 'absmin']]
stats_row_avg = stats_pickle2.groupby('row').mean()
stats_row_avg['sensor_loc'] = 'avg'
# stats_row_avg = stats_row_avg.set_index('sensor_loc',append=True)
# #stats_total = pd.concat([stats_pickle2, stats_row_avg])
# stats_total = stats_pickle2.merge(stats_row_avg, left_index = True, right_index=True, how="inner") 

# stats_row_avg.to_csv('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats.csv')
stats_row_avg.to_csv('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats_no_adj.csv')

#%%
# new data frame with split value columns
new = median_pickle["sensor"].str.split("R|_|_", n = 3, expand = True)
 
# making separate first name column from new data frame
median_pickle["row"]= new[1]
 
# making separate last name column from new data frame
median_pickle["sensor_loc"]= new[2]

median_pickle = median_pickle.rename(columns={'trough_angle_dev':'tracking_error'})

median_row_avg = median_pickle.groupby(['time','season','row']).mean(numeric_only=True)

# median_row_avg.to_csv('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_median_day.csv')
median_row_avg.to_csv('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_median_day_no_adj.csv')

#%%
intercept_pickle2 = intercept_pickle['DO']['intercept_factor']
incercept_row_avg = intercept_pickle2.groupby(['row','stat']).mean()
incercept_row_avg.to_csv('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/intercept_factor_stats.csv')
