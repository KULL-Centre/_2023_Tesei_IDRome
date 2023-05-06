import numpy as np
import pandas as pd
from localcider.sequenceParameters import SequenceParameters
import glob
import difflib
from sklearn import svm
from argparse import ArgumentParser
import itertools
from omadb import Client
from joblib import dump, load

parser = ArgumentParser()
parser.add_argument('--batch',nargs='?',const='', type=int)
args = parser.parse_args()
batch = args.batch

def predict_nu(i,residues,filename,model,features):
    omadb_client = Client()
    df = pd.read_csv(filename,
                    header=None,delimiter='>',index_col=1,names=['fasta','protein'])
    print(filename)
    if df.index.isnull().sum() > 1 and len(df.index[0].split('to')) > 1:
        df.loc[~df.index.isnull()] = df.loc[df.index.isnull()].values
        df = df.loc[~df.index.isnull()]
        df['fasta'] = df.loc[~df.index.isnull()].fasta.apply(lambda x: x.replace('-',''))
        df['N'] = df.fasta.apply(len)
        df.reset_index(inplace=True)
        human_id = filename.split('/')[1].split('.')[0]
        df['HUMAN'] = human_id
        df[['protein','first','to','last']] = df.protein.str.split(' ', 3, expand=True)
        df.set_index('protein',inplace=True,drop=False)
        df.loc[df.index.str.startswith(human_id),'protein'] = human_id
        df.loc[df.index.str.startswith(human_id),'first'] = human_id.split('_')[1].split('to')[0]
        df.loc[df.index.str.startswith(human_id),'last'] = human_id.split('to')[1]
        df.loc[~df.index.str.startswith(human_id),'protein'] = df.loc[~df.index.str.startswith(human_id),'protein'].apply(
                    lambda x: x + '_{:s}to{:s}'.format(df.loc[x,'first'],df.loc[x,'last']))
        df = df.drop('to',axis=1).set_index('protein')
        for homologue in df.index:
            print(homologue)
            fasta_homologue = df.loc[homologue].fasta
            N = df.loc[homologue].N
            if (N >= 30) and np.all([res in residues.index for res in fasta_homologue]):
                df.loc[homologue,'mean_lambda'] = np.mean(residues.loc[list(fasta_homologue)].lambdas)
                df.loc[homologue,'faro'] = sum([fasta_homologue.count(a) for a in ['W','Y','F']])/N
                df.loc[homologue,'fR'] = fasta_homologue.count('R')/N
                df.loc[homologue,'fK'] = fasta_homologue.count('K')/N
                df.loc[homologue,'fE'] = fasta_homologue.count('E')/N
                df.loc[homologue,'fD'] = fasta_homologue.count('D')/N

                pairs = np.array(list(itertools.combinations(fasta_homologue,2)))
                pairs_indices = np.array(list(itertools.combinations(range(N),2)))
                # calculate sequence separations
                ij_dist = np.diff(pairs_indices,axis=1).flatten().astype(float)
                # calculate lambda sums
                ll = residues.lambdas.loc[pairs[:,0]].values+residues.lambdas.loc[pairs[:,1]].values
                # calculate SHD
                beta = -1
                df.loc[homologue,'shd'] = np.sum(ll*np.power(np.abs(ij_dist),beta))/N

                N_FL = 0
                omaid = ''
                canonicalid = ''
                try:
                    omadb_dict = omadb_client.proteins[homologue.split('_')[0]]
                    N_FL = omadb_dict['sequence_length']
                    omaid = omadb_dict['omaid']
                    canonicalid = omadb_dict['canonicalid']
                except:
                    omadb_dict = omadb_client.proteins.search(fasta_homologue)
                    if len(omadb_dict['targets']) > 1:
                        seq_lengths = []
                        for targets_k in omadb_dict['targets']:
                            seq_lengths.append(targets_k['sequence_length'])
                            if targets_k['omaid'] == homologue.split('_')[0]:
                                N_FL = targets_k['sequence_length']
                                omaid = targets_k['omaid']
                                canonicalid = targets_k['canonicalid']
                        if N_FL == 0 and omaid == '' and canonicalid == '':
                            N_FL = omadb_dict['targets'][np.argmax(seq_lengths)]['sequence_length']
                            omaid = omadb_dict['targets'][np.argmax(seq_lengths)]['omaid']
                            canonicalid = omadb_dict['targets'][np.argmax(seq_lengths)]['canonicalid']

                df.loc[homologue,'N_FL'] = N_FL
                df.loc[homologue,'omaid'] = omaid
                df.loc[homologue,'canonicalid'] = canonicalid

                fasta = list(fasta_homologue).copy()
                fasta_kappa = fasta.copy()
                if int(df.loc[homologue,'first']) < 2:
                    residues.loc['X'] = residues.loc[fasta_homologue[0]]
                    residues.loc['X','q'] = residues.loc[fasta_homologue[0],'q'] + 1.
                    fasta[0] = 'X'
                    if residues.loc['X','q'] > 0:
                        fasta_kappa[0] = 'K'
                    else:
                        fasta_kappa[0] = 'A'
                if int(df.loc[homologue,'last']) >= N_FL-1+int(df.loc[homologue,'first']):
                    residues.loc['Z'] = residues.loc[fasta_homologue[-1]]
                    residues.loc['Z','q'] = residues.loc[fasta_homologue[-1],'q'] - 1.
                    fasta[-1] = 'Z'
                    if residues.loc['Z','q'] < 0:
                        fasta_kappa[-1] = 'D'
                    else:
                        fasta_kappa[-1] = 'A'

                pairs = np.array(list(itertools.combinations(fasta,2)))
                # calculate charge products
                qq = residues.q.loc[pairs[:,0]].values*residues.q.loc[pairs[:,1]].values
                # calculate SCD
                df.loc[homologue,'scd'] = np.sum(qq*np.sqrt(ij_dist))/N

                SeqOb = SequenceParameters(''.join(fasta_kappa))
                kappa = SeqOb.get_kappa()
                df.loc[homologue,'kappa'] = 0 if kappa==-1 else kappa
                df.loc[homologue,'fcr'] = residues.q.loc[fasta].abs().mean()
                df.loc[homologue,'ncpr'] = residues.q.loc[fasta].mean()
                df.loc[homologue,'nu'] = model.predict(df.loc[homologue,features].values.reshape(1, -1))
        return df

residues = pd.read_csv('../../md_simulations/data/residues.csv').set_index('one',drop=False)

features = ['scd','shd','kappa','fcr','mean_lambda']
model = load('../../svr_model/svr_model.joblib')

# The set of homologs was downloaded from https://doi.org/10.5281/zenodo.6311384
# https://zenodo.org/record/6311384/files/human_idr_homologues.zip
homology_folder = 'human_idr_homologues'

dfs = []
for i,filename in enumerate(glob.glob(homology_folder+'/*')[batch:batch+1000]):
    dfs.append(predict_nu(i,residues,filename,model,features))

pd.concat(dfs).to_csv(f'nu_svm_{batch:d}.csv')

