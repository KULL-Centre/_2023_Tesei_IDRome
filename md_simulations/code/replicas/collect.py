import subprocess
import os
import pandas as pd
import numpy as np
import mdtraj as md
import time
from jinja2 import Template
from statsmodels.tsa.stattools import acf

sequences = pd.read_csv('IDRome_DB.csv',index_col=0)

for name in ['Q6IMN6_317_998','Q5VUA4_1_884','P46100_289_1546','Q15648_525_1581','Q13342_131_587','Q9UQL6_176_680']:
    print(name,sequences.loc[name].fasta)
    nres = len(sequences.loc[name].fasta)
    nu = []
    R0 = []
    ete2_Rg2 = []
    emap = []
    stdev = []
    rgarray = []
    Delta_array = []
    S_array = []
    acf_1 = []
    acf_2 = []
    acf_3 = []
    Delta = []
    S = []
    SPR = []
    for replica in range(5):
        analysis = pd.read_csv(name+f'/{replica:d}/analysis.csv',index_col=0)
        emap.append(pd.read_csv(name+f'/{replica:d}/emap.csv',index_col=0))
        stdev.append(pd.read_csv(name+f'/{replica:d}/stdev.csv',index_col=0))
        Delta.append(analysis.loc['Delta','value'])
        S.append(analysis.loc['S','value'])
        SPR.append(analysis.loc['SPR','value'])
        nu.append(analysis.loc['nu','value'])
        R0.append(analysis.loc['R0','value'])
        ete2_Rg2.append(analysis.loc['ete2_Rg2','value'])
        rg = np.load(name+f'/{replica:d}/rg.npy')
        rgarray.append(rg)
        Delta_array.append(np.load(name+f'/{replica:d}/Delta.npy'))
        S_array.append(np.load(name+f'/{replica:d}/S.npy'))
        acf_rg = acf(rg,nlags=5,fft=True)
        acf_1.append(acf_rg[1])
        acf_2.append(acf_rg[2])
        acf_3.append(acf_rg[3])
    pd.DataFrame(index=range(nres),columns=range(nres),data=np.mean(emap,axis=0)).to_csv(name+'_emap.csv.gz')
    pd.DataFrame(index=range(nres),columns=range(nres),data=np.mean(stdev,axis=0)).to_csv(name+'_stdev.csv.gz')
    sequences.loc[name,'acf_rg_1'] = np.mean(acf_1)
    sequences.loc[name,'acf_rg_1_err'] = np.std(acf_1)
    sequences.loc[name,'acf_rg_2'] = np.mean(acf_2)
    sequences.loc[name,'acf_rg_2_err'] = np.std(acf_2)
    sequences.loc[name,'acf_rg_3'] = np.mean(acf_3)
    sequences.loc[name,'acf_rg_3_err'] = np.std(acf_3)
    sequences.loc[name,'Delta'] = np.mean(Delta)
    sequences.loc[name,'Delta_err'] = np.std(Delta)
    sequences.loc[name,'S'] = np.mean(S)
    sequences.loc[name,'S_err'] = np.std(S)
    sequences.loc[name,'SPR'] = np.mean(SPR)
    sequences.loc[name,'SPR_err'] = np.std(SPR)
    sequences.loc[name,'nu_replicas'] = ', '.join([str(x) for x in nu])
    sequences.loc[name,'nu'] = np.mean(nu)
    sequences.loc[name,'nu_err'] = np.std(nu)
    sequences.loc[name,'R0_replicas'] = ', '.join([str(x) for x in R0])
    sequences.loc[name,'R0'] = np.mean(R0)
    sequences.loc[name,'R0_err'] = np.std(R0)
    sequences.loc[name,'ete2_Rg2_replicas'] = ', '.join([str(x) for x in ete2_Rg2])
    sequences.loc[name,'ete2_Rg2'] = np.mean(ete2_Rg2)
    sequences.loc[name,'ete2_Rg2_err'] = np.std(ete2_Rg2)
    np.save(name+'_rgarray.npy',np.array(rgarray).flatten())
    np.save(name+'_Delta_array.npy',np.array(Delta_array).flatten())
    np.save(name+'_S_array.npy',np.array(S_array).flatten())
    print(name)
    time.sleep(.6)

sequences.to_csv('replicas_data.csv')
