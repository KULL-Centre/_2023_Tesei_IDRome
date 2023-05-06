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
#SBATCH -t 72:00:00
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /groups/sbinlab/giulio/.bashrc

conda activate idrome

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --seq_name {{name}} --path {{path}}""")

sequences = pd.read_csv('../../data/replicas_data.csv',index_col=0)

for name in sequences.index:
    print(name)
    if not os.path.isdir(name):
        os.mkdir(name)
    for replica in range(5):
        if not os.path.isdir(name+'/{:d}'.format(replica)):
            os.mkdir(name+'/{:d}'.format(replica))
        with open('{:s}_{:d}.sh'.format(name,replica), 'w') as submit:
            submit.write(submission.render(name=name,replica=replica,path=name+'/{:d}'.format(replica)))
        subprocess.run(['sbatch','{:s}_{:d}.sh'.format(name,replica)])
        print(name)
        time.sleep(.6)
