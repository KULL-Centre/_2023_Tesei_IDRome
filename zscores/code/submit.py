import pandas as pd
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH --job-name={{batch}}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=2GB
#SBATCH -t 100:00:00
#SBATCH -o {{batch}}_out
#SBATCH -e {{batch}}_err

source /groups/sbinlab/giulio/.bashrc

conda activate nardini

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./analyse.py --batch {{batch}}
""")


df_idrome = pd.read_csv('IDRome_DB.csv',index_col=0)

for batch in range(0,df_idrome.shape[0],50):
    with open(f'{batch:d}.sh', 'w') as submit:
        submit.write(submission.render(batch=batch))
    subprocess.run(['sbatch',f'{batch:d}.sh'])
