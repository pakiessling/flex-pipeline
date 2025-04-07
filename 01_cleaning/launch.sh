#!/usr/local_rwth/bin/zsh

nohup snakemake --conda-frontend conda --profile slurm --jobscript  jobscript.sh --config samples="samples.csv" > pipeline.log &
