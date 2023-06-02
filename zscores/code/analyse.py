import pandas as pd
from argparse import ArgumentParser
from nardini.constants import TYPEALL
from nardini.utils import read_sequences_from_string_list
from nardini.score_and_plot import calculate_zscore
import time

parser = ArgumentParser()
parser.add_argument('--batch',nargs='?',const='', type=int)
args = parser.parse_args()
batch = args.batch

def calc_z_scores(df_idrome,seq_name,r):
    fasta = list(df_idrome.loc[seq_name].fasta)
    fasta_kappa = fasta.copy()

    # calculate properties that do not depend on charges
    fasta_sequence = read_sequences_from_string_list([''.join(fasta)],seq_name)
    data_to_export = calculate_zscore(fasta_sequence, TYPEALL, int(5e5), 1684533501)
    df = pd.DataFrame(data_to_export).T.rename(
        {0:'seq',1:'seq_scr',2:'index',3:'z_mat',4:'z_mat_scr'},axis=1)
    df.rename(index={df.index.values[0]: seq_name},inplace=True)
    df['seq_name'] = seq_name
    df['termini'] = False
    df['delta_+-'] = df.z_mat.values[0][2,3]
    df['omega_h'] = df.z_mat.values[0].diagonal()[1]
    df['omega_+'] = df.z_mat.values[0].diagonal()[2]
    df['omega_-'] = df.z_mat.values[0].diagonal()[3]
    df['omega_pi'] = df.z_mat.values[0].diagonal()[4]
    df.drop('index',axis=1,inplace=True)
    # fix charges
    if df_idrome.loc[seq_name,'first'] == 1:
        r.loc['X'] = r.loc[fasta[0]]
        r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
        fasta[0] = 'X'
        if r.loc['X','q'] > 0:
            fasta_kappa[0] = 'K'
        else:
            fasta_kappa[0] = 'A'
    if df_idrome.loc[seq_name,'last'] == df_idrome.loc[seq_name,'N_FL']:
        r.loc['Z'] = r.loc[fasta[-1]]
        r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
        fasta[-1] = 'Z'
        if r.loc['Z','q'] < 0:
            fasta_kappa[-1] = 'D'
        else:
            fasta_kappa[-1] = 'A'

    if fasta_kappa != fasta:
        # calculate properties that depend on charges
        fasta_sequence = read_sequences_from_string_list([''.join(fasta_kappa)],seq_name+'_termini')
        data_to_export = calculate_zscore(fasta_sequence, TYPEALL, int(5e5), 1684533501)
        df_termini = pd.DataFrame(data_to_export).T.rename(
            {0:'seq',1:'seq_scr',2:'index',3:'z_mat',4:'z_mat_scr'},axis=1)
        df_termini.rename(index={df_termini.index.values[0]: seq_name+'_termini'},inplace=True)
        df_termini['seq_name'] = seq_name
        df_termini['termini'] = True
        df_termini['delta_+-'] = df_termini.z_mat.values[0][2,3]
        df_termini['omega_h'] = df_termini.z_mat.values[0].diagonal()[1]
        df_termini['omega_+'] = df_termini.z_mat.values[0].diagonal()[2]
        df_termini['omega_-'] = df_termini.z_mat.values[0].diagonal()[3]
        df_termini['omega_pi'] = df_termini.z_mat.values[0].diagonal()[4]
        df_termini.drop('index',axis=1,inplace=True)
        df = pd.concat([df,df_termini])
    return df

t0 = time.time()

df_idrome = pd.read_csv('IDRome_DB.csv',index_col=0)
r = pd.read_csv('residues.csv').set_index('one',drop=False)

dfs = []
for seq_name in df_idrome.iloc[batch:batch+50].index:
    dfs.append(calc_z_scores(df_idrome,seq_name,r))
    print(seq_name)

pd.concat(dfs).to_csv(f'zscores_{batch:d}.csv.gz')

print('Time: {:.3f} h'.format((time.time()-t0)/3600))
