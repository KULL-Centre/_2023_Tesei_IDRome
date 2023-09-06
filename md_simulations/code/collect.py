import subprocess
import os
import pandas as pd
import numpy as np
import time
import mdtraj as md

sequences = pd.read_csv('sequences.csv',index_col=0)

df_idrome = pd.read_csv('conf_prop_long_idrs.csv',index_col=0)
print(df_idrome.shape,df_idrome.dropna().shape)
df_idrome = df_idrome[df_idrome.index.isin(sequences.index)]

residues = pd.read_csv('residues.csv').set_index('one',drop=False)

f = lambda x,R0,v : R0*np.power(x,v)

for j,name in enumerate(sequences.iloc[:].index):
    first_last = '{:g}_{:g}'.format(sequences.loc[name,'first'],sequences.loc[name,'last'])
    title1 = name[:2]
    title2 = name[2:4]
    title3 = name[4:6]
    title4 = name[6:10]
    if len(sequences.loc[name,'uniprot']) == 6:
        path = os.path.join(title1,title2,title3,first_last)
    else:
        path = os.path.join(title1,title2,title3,title4,first_last)
    if os.path.isfile(path+'/analysis.csv'):
        df_analysis = pd.read_csv(path+'/analysis.csv',index_col=0)
        for prop in ['nu','S','Delta','ete2_Rg2','Rg','ete','rh','R0','SPR']:
            df_idrome.loc[name,prop] = df_analysis.loc[prop,'value']
            df_idrome.loc[name,prop+'_err'] = df_analysis.loc[prop,'error']

print((~df_idrome.nu.isnull()).sum())

df_idrome.to_csv('conf_prop.csv.gz')

print(sequences.shape[0],df_idrome.dropna().shape[0],sequences.loc[df_idrome.index].nres_seg.max())
