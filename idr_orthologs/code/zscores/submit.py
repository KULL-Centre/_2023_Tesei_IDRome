import pandas as pd
import subprocess
import time
import os
from jinja2 import Template

submission = Template("""#!/bin/sh
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N {{batch}}
#PBS -m n
#PBS -l nodes=1:ppn=1:thinnode
#PBS -l mem=1gb
#PBS -l walltime=120:00:00
#PBS -e {{batch}}_err
#PBS -o {{batch}}_out

source /home/people/giutes/.bashrc
conda activate idrome

cd $PBS_O_WORKDIR

python ./calc_seq_prop.py --batch {{batch}}
""")

df_o = pd.read_csv(f'idr_orthologs.csv.gz',index_col=0,dtype='object',comment='#')
df_o = df_o.reset_index()
df_o = df_o.rename({'ortholog_full_len':'N_FL','idr_full_len':'human_N_FL',
                              'idr':'human','ortholog_seq':'fasta','idr_ortholog':'seq_name'},axis=1)
df_o = df_o.set_index('seq_name')
df_o['first'] = [int(i.split('_')[1]) for i in df_o.index]
df_o['last'] = [int(i.split('_')[2]) for i in df_o.index]
df_o = df_o[~df_o.index.duplicated(keep='first')]

for batch in range(0,df_o.shape[0],20):
    if not os.path.isfile(f'seq_prop_{batch:d}.csv.gz'):
        with open(f'{batch:d}.sh', 'w') as submit:
            submit.write(submission.render(batch=batch))
        subprocess.run(['qsub',f'{batch:d}.sh'])
    time.sleep(.6)
