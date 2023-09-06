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
#SBATCH --mem=30GB
#SBATCH -t 200:00:00
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --seq_name {{name}} --path {{path}}""")

sequences = pd.read_csv('rg_test_set.csv',index_col=0)

# Download blocking code
if not os.path.exists('BLOCKING'):
    subprocess.check_call(['git','clone','https://github.com/fpesceKU/BLOCKING'])

for name in sequences.index:
    print(name)
    if not os.path.isdir(name):
        os.mkdir(name)
        with open(f'{name:s}.sh', 'w') as submit:
            submit.write(submission.render(name=name,path=name))
        subprocess.run(['sbatch',f'{name:s}.sh'])
        print(name)
        time.sleep(.6)
