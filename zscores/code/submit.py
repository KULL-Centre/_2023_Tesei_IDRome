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
#PBS -l walltime=72:00:00
#PBS -e {{batch}}_err
#PBS -o {{batch}}_out

source /home/people/giutes/.bashrc
conda activate idrome

cd $PBS_O_WORKDIR

python ./calc_seq_prop.py --batch {{batch}}
""")

df_idrome = pd.read_csv('IDRome_DB.csv',index_col=0)

for batch in range(0,df_idrome.shape[0],20):
    if not os.path.isfile(f'seq_prop_{batch:d}.csv.gz'):
        print(batch)
        with open(f'{batch:d}.sh', 'w') as submit:
            submit.write(submission.render(batch=batch))
        subprocess.run(['qsub',f'{batch:d}.sh'])
    time.sleep(.6)
