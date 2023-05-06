import os
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/bash
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N {{batch}}
#PBS -m n
#PBS -l nodes=1:ppn=1:thinnode
#PBS -l mem=1gb
#PBS -l walltime=10:00:00
#PBS -e {{batch}}.err
#PBS -o {{batch}}.out

source /home/people/giutes/.bashrc
conda activate idrome

cd $PBS_O_WORKDIR

python ./analyse.py --batch {{batch}}
""")


for batch in range(0,15996,1000):
    with open(f'{batch:d}.sh', 'w') as submit:
        submit.write(submission.render(batch=batch))
    subprocess.run(['qsub',f'{batch:d}.sh'])
