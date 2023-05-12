import subprocess
import os
import pandas as pd
import numpy as np
import mdtraj as md
import time
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH --job-name={{name}}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=2GB
#SBATCH -t 10:00:00
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /groups/sbinlab/giulio/.bashrc

conda activate idrome

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --seq_name {{name}} --path {{path}}""")

sequences = pd.read_csv('../data/idr_all.csv.gz',header=0,sep=';')
sequences.sort_values('uniprot',inplace=True)
sequences['seq_name'] = sequences.uniprot+'_'+sequences['first'].apply(lambda x : '{:g}'.format(x))+'_'+sequences['last'].apply(lambda x : '{:g}'.format(x))
sequences.set_index('seq_name',inplace=True)

# Download blocking code
if not os.path.exists('BLOCKING'):
    subprocess.check_call(['git','clone','https://github.com/fpesceKU/BLOCKING'])

print(sequences.index.size)
# Making the directories, sorted by their uniprot IDs split into twos, i.e. P57678_180_217 will be saved to P5/76/78/180_217
for i,name in enumerate(sequences.index[:]):
    first_last = '{:g}_{:g}'.format(sequences.loc[name,'first'],sequences.loc[name,'last'])
    title1 = name[:2]
    title2 = name[2:4]
    title3 = name[4:6]
    title4 = name[6:10]
    if len(sequences.loc[name,'uniprot']) == 6:
    	path = os.path.join(title1,title2,title3,first_last)
    else:
    	path = os.path.join(title1,title2,title3,title4,first_last)
    print(path)
    if not os.path.exists(path):
        os.makedirs(path)
        with open('{:s}.sh'.format(name), 'w') as submit:
            submit.write(submission.render(name=name,path=path))
        subprocess.run(['sbatch','{:s}.sh'.format(name)])
        print(i,path)
    elif np.load(path+'/rg.npy').size != 1000:
        with open('{:s}.sh'.format(name), 'w') as submit:
            submit.write(submission.render(name=name,path=path))
        subprocess.run(['sbatch','{:s}.sh'.format(name)])
        print(i,path)
    time.sleep(.6)
