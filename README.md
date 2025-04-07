# 10X FLEX processing pipeline

Change the Conda path in 00_cleaning/jobscript.sh and 01_integration/jobscript.sh
*TODO fix this with executor plugin refactor* 

## Step 1 cleaning
- Fill in the paths to the two Cellranger .h5 outputs for every sample in 00_cleaning/samples.csv
- Execute 00_cleaning/launch.sh (needs Snakemake in the PATH and SLURM)

## Step 2 integration
- When you are satisfied with the cleaning add the path to a folder containing .h5ad files to 01_integration/launch.sh
- Execute 01_integration/launch.sh (needs Snakemake in the PATH and SLURM)