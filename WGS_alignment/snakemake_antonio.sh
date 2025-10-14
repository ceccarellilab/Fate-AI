#Snakemake Slurm 9.12.0
mamba create -c conda-forge -c bioconda -c nodefaults -n snakemake_last snakemake snakemake-executor-plugin-slurm
conda activate snakemake_last
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake

