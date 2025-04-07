#!/usr/local_rwth/bin/zsh



if [ -d "$HOME/miniconda3" ] && [ -x "$HOME/miniconda3" ]; then
    export CONDA_ROOT="$HOME/miniconda3"
else
    export CONDA_ROOT="/work/mz637064/miniconda3/" # path to your conda or miniconda3
fi


. $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"

source activate /hpcwork/p0020567/enviroments/sopa2 # path to your Snakemake env

nohup snakemake --executor slurm --conda-frontend conda --config input="00_cleaning/output"   > pipeline.log &
