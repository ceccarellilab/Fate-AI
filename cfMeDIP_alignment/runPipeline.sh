#snakemake -s cfMedip_split.Snakefile --use-conda --jobs 17 --resources total_cpus=220 total_mem_mb=312000 --latency-wait 320 --keep-going --cluster-config cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}"


CONFIG_FILE=""
FIXED_CONDA_PATH="/home2/adefalco/Fate-AI/output_folder/.snakemake/conda/"
SAMPLE_IN_PARALLEL="4"

## ALIGNMENT
snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s cfMedip_split.Snakefile \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
