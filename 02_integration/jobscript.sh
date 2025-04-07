#!/usr/local_rwth/bin/zsh

if [ -d "$HOME/miniconda3" ] && [ -x "$HOME/miniconda3" ]; then
    export CONDA_ROOT="$HOME/miniconda3"
else
    export CONDA_ROOT="/work/mz637064/miniconda3/" # path to your conda or miniconda3
fi


. $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"

{exec_job}
