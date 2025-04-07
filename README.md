# 10X FLEX processing pipeline

Change the Conda path ```in 01_cleaning/jobscript.sh``` and ```02_integration/jobscript.sh```


*TODO fix this with executor plugin refactor* 

## Step 1 cleaning
- Fill in the paths to the two Cellranger .h5 outputs for every sample in ```01_cleaning/samples.csv```

  i.e *sample_filtered_feature_bc_matrix.h5 and sample_raw_feature_bc_matrix.h5*
- Execute ```01_cleaning/launch.sh``` (needs Snakemake in the PATH and SLURM)

## Step 2 integration
- When you are satisfied with the cleaning add the path to a folder containing .h5ad files to ```02_integration/launch.sh```
- Execute ```02_integration/launch.sh``` (needs Snakemake in the PATH and SLURM)
