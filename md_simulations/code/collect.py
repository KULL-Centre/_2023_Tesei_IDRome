import subprocess
import os
import pandas as pd
import numpy as np
import time

sequences = pd.read_csv('../data/idr_all.csv.gz',header=0,sep=';')
sequences.sort_values('uniprot',inplace=True)
sequences['seq_name'] = sequences.uniprot+'_'+sequences['first'].apply(lambda x : '{:g}'.format(x))+'_'+sequences['last'].apply(lambda x : '{:g}'.format(x))
sequences.set_index('seq_name',inplace=True)

df_idrome = pd.DataFrame(index=sequences.index)

for i,name in enumerate(sequences.iloc[:].index):
    first_last = '{:g}_{:g}'.format(sequences.loc[name,'first'],sequences.loc[name,'last'])
    title1 = name[:2]
    title2 = name[2:4]
    title3 = name[4:6]
    title4 = name[6:10]
    if len(sequences.loc[name,'uniprot']) == 6:
    	path = os.path.join(title1,title2,title3,first_last)
    else:
    	path = os.path.join(title1,title2,title3,title4,first_last)
    df_analysis = pd.read_csv(path+'/analysis.csv',index_col=0)
    for prop in ['nu','S','Delta','ete2_Rg2','Rg','ete','rh','R0']:
        df_idrome.loc[name,prop] = df_analysis.loc[prop,'value']
        df_idrome.loc[name,prop+'_err'] = df_analysis.loc[prop,'error']
df_idrome.to_csv('../data/conf_prop.csv')
