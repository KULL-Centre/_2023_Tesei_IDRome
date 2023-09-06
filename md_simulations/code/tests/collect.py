import subprocess
import os
import pandas as pd
import numpy as np
import time

sequences = pd.read_csv('rg_test_set.csv',index_col=0)
for name in sequences.index:
    df_analysis = pd.read_csv(name+'/analysis.csv',index_col=0)
    for prop in ['nu','S','Delta','ete2_Rg2','Rg','ete','rh','R0']:
        sequences.loc[name,prop] = df_analysis.loc[prop,'value']
        sequences.loc[name,prop+'_err'] = df_analysis.loc[prop,'error']
sequences.to_csv('conf_rg_test_set.csv')
