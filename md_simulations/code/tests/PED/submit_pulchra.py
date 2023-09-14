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
#SBATCH --exclusive
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=20GB
#SBATCH -t 72:00:00
#SBATCH -o all_atom/{{name}}_out
#SBATCH -e all_atom/{{name}}_err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./pulchra.py --name {{name}} --num_cpus 32 --pulchra /sbinlab/giulio/pulchra_306/pulchra""")

fret_IDRs = ['NLS','NUS','IBB','NUL','Nup49','Sic1']
sequences = pd.read_csv('PED_sequences.csv',index_col=0)

for name in fret_IDRs:
    with open('{:s}_pulchra.sh'.format(name), 'w') as submit:
        submit.write(submission.render(name=name))
    subprocess.run(['sbatch','{:s}_pulchra.sh'.format(name)])
    print(name)
    time.sleep(.6)
