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
#SBATCH --partition=qgpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=140GB
#SBATCH -t 48:00:00
#SBATCH -e {{path}}/err
#SBATCH -o {{path}}/out

source /home/gitesei/.bashrc
conda activate idrome
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python ./simulate.py --seq_name {{name}} --path {{path}}""")

df_prefilter = pd.read_csv('idr.csv.gz',header=0,sep=';')
df_prefilter['seq_name'] = df_prefilter.uniprot+'_'+df_prefilter['first'].apply(lambda x :
                        '{:g}'.format(x))+'_'+df_prefilter['last'].apply(lambda x : '{:g}'.format(x))
df_prefilter.set_index('seq_name',inplace=True)
sequences = sequences.sort_values('nres_seg',ascending=False)
sequences.to_csv('sequences.csv')

# Download blocking code
if not os.path.exists('BLOCKING'):
    subprocess.check_call(['git','clone','https://github.com/fpesceKU/BLOCKING'])

print(sequences.index.size,df_long_idrs.index.size)
#Making the directories, sorted by their uniprot IDs split into twos, i.e. P57678_180_217 will be saved to P5/76/78/180_217
for i,name in enumerate(sequences.index[:]):
    if name not in df_long_idrs.index:
        first_last = '{:g}_{:g}'.format(sequences.loc[name,'first'],sequences.loc[name,'last'])
        title1 = name[:2]
        title2 = name[2:4]
        title3 = name[4:6]
        title4 = name[6:10]
        if len(sequences.loc[name,'uniprot']) == 6:
            path = os.path.join(title1,title2,title3,first_last)
        else:
            path = os.path.join(title1,title2,title3,title4,first_last)
        #if not os.path.exists(path):
        if not os.path.isfile(path+'/traj.xtc'):
            if not os.path.exists(path):
                os.makedirs(path)
            print(path)
            with open('{:s}.sh'.format(name), 'w') as submit:
                submit.write(submission.render(name=name,path=path))
            subprocess.run(['sbatch','{:s}.sh'.format(name)])
