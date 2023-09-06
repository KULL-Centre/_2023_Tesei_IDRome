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
#SBATCH -t 5:00:00
#SBATCH -e {{path}}/err
#SBATCH -o {{path}}/out

source /home/gitesei/.bashrc
conda activate idrome
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python ./simulate.py --seq_name {{name}} --path {{path}}""")

sequences = pd.read_csv('conf_buffering_seq.csv',index_col=0)

# Download blocking code
if not os.path.exists('BLOCKING'):
    subprocess.check_call(['git','clone','https://github.com/fpesceKU/BLOCKING'])

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
