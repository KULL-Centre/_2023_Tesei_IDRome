import pandas as pd
from argparse import ArgumentParser
from nardini.constants import TYPEALL
from nardini.utils import read_sequences_from_string_list
from nardini.score_and_plot import calculate_zscore
from localcider.sequenceParameters import SequenceParameters
import time
import numpy as np
import itertools
from joblib import dump, load

parser = ArgumentParser()
parser.add_argument('--batch',nargs='?',const='', type=int)
args = parser.parse_args()
batch = args.batch

def calc_seq_prop(df_o,seq_name,r,model_nu,model_spr):
    columns_z_mat = np.array(list(itertools.product(['mu','h',
        '+','-','pi','A','P','G'],repeat=2)))
    columns_z_mat = np.apply_along_axis('_'.join, axis=1,
            arr=columns_z_mat).reshape(8,8)
    columns_triu = columns_z_mat[np.triu_indices(8)]
    # input features for SVR models
    features_nu = ['scd','shd','kappa','fcr','mean_lambda']
    features_spr = ['scd','shd','mean_lambda']

    fasta = list(df_o.loc[seq_name].fasta)
    fasta_kappa = fasta.copy()
    N = len(fasta)

    # calculate properties that do not depend on charges
    fK = sum([fasta.count(a) for a in ['K']])/N
    fR = sum([fasta.count(a) for a in ['R']])/N
    fE = sum([fasta.count(a) for a in ['E']])/N
    fD = sum([fasta.count(a) for a in ['D']])/N
    mean_lambda = np.mean(r.loc[fasta].lambdas)
    pairs = np.array(list(itertools.combinations(fasta,2)))
    pairs_indices = np.array(list(itertools.combinations(range(N),2)))
    # calculate sequence separations
    ij_dist = np.diff(pairs_indices,axis=1).flatten().astype(float)
    # calculate lambda sums
    ll = r.lambdas.loc[pairs[:,0]].values+r.lambdas.loc[pairs[:,1]].values
    # calculate SHD
    beta = -1
    shd = np.sum(ll*np.power(np.abs(ij_dist),beta))/N

    fasta_sequence = read_sequences_from_string_list([''.join(fasta)],seq_name)
    data_to_export = calculate_zscore(fasta_sequence, TYPEALL, int(1e5), 1684533501)
    df = pd.DataFrame(data_to_export).T.rename(
        {0:'seq',1:'seq_scr',2:'index',3:'z_mat',4:'z_mat_scr'},axis=1)
    df.rename(index={df.index.values[0]: seq_name},inplace=True)
    df['seq_name'] = seq_name
    df['termini'] = False
    df[columns_triu] = pd.Series(data=df.z_mat.values[0].flatten(),
            index=columns_z_mat.flatten()).loc[columns_triu]
    df.drop('index',axis=1,inplace=True)

    # fix charges
    if df_o.loc[seq_name,'first'] == 1:
        r.loc['X'] = r.loc[fasta[0]]
        r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
        if r.loc['X','q'] > 0:
            if r.loc[fasta[0],'q'] == 0:
                fasta_kappa[0] = 'K'
        else:
            fasta_kappa[0] = 'A'
        fasta[0] = 'X'
    if df_o.loc[seq_name,'last'] == df_o.loc[seq_name,'N_FL']:
        r.loc['Z'] = r.loc[fasta[-1]]
        r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
        if r.loc['Z','q'] < 0:
            if r.loc[fasta[-1],'q'] == 0:
                fasta_kappa[-1] = 'D'
        else:
            fasta_kappa[-1] = 'A'
        fasta[-1] = 'Z'

    # calculate properties that depend on charges
    pairs = np.array(list(itertools.combinations(fasta,2)))
    # calculate charge products
    qq = r.q.loc[pairs[:,0]].values*r.q.loc[pairs[:,1]].values
    # calculate SCD
    scd = np.sum(qq*np.sqrt(ij_dist))/N
    SeqOb = SequenceParameters(''.join(fasta_kappa))
    kappa = SeqOb.get_kappa()
    kappa = 0 if kappa==-1 else kappa
    fcr = r.q.loc[list(fasta)].abs().mean()
    ncpr = r.q.loc[list(fasta)].mean()

    df[['fK','fR','fE','fD','mean_lambda','shd',
        'scd','kappa','fcr','ncpr']] = [fK,fR,fE,fD,mean_lambda,shd,scd,kappa,fcr,ncpr]
    nu_svr = float(model_nu.predict(df[features_nu].values.reshape(1, -1)))
    spr_svr = float(model_spr.predict(df[features_spr].values.reshape(1, -1)))
    df['nu_svr'] = nu_svr
    df['SPR_svr'] = spr_svr

    if fasta_kappa != list(df_o.loc[seq_name].fasta):
        fasta_sequence = read_sequences_from_string_list(
                [''.join(fasta_kappa)],seq_name+'_termini')
        data_to_export = calculate_zscore(fasta_sequence, TYPEALL, int(1e5), 1684533501)
        df_termini = pd.DataFrame(data_to_export).T.rename(
            {0:'seq',1:'seq_scr',2:'index',3:'z_mat',4:'z_mat_scr'},axis=1)
        df_termini.rename(index={df_termini.index.values[0]: seq_name+'_termini'},
                inplace=True)
        df_termini['seq_name'] = seq_name
        df_termini['termini'] = True
        df_termini[columns_triu] = pd.Series(data=df_termini.z_mat.values[0].flatten(),
                index=columns_z_mat.flatten()).loc[columns_triu]
        df_termini.drop('index',axis=1,inplace=True)
        df_termini[['fK','fR','fE','fD','mean_lambda','shd',
            'scd','kappa','fcr','ncpr','nu_svr','SPR_svr']] = [fK,fR,fE,fD,mean_lambda,shd,
            scd,kappa,fcr,ncpr,nu_svr,spr_svr]
        df = pd.concat([df,df_termini])
    return df

t0 = time.time()

df_idrome = pd.read_csv('IDRome_DB.csv',index_col=0)

r = pd.read_csv('residues.csv').set_index('one',drop=False)

model_nu = load('svr_model_nu.joblib')
model_spr = load('svr_model_SPR.joblib')

print('Starting the loop')

dfs = []
for seq_name in df_idrome.iloc[batch:batch+20].index:
    dfs.append(calc_seq_prop(df_idrome,seq_name,r,model_nu,model_spr))
    print(seq_name)

pd.concat(dfs).to_csv(f'seq_prop_{batch:d}.csv.gz')

print('Time: {:.3f} h'.format((time.time()-t0)/3600))


