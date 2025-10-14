CONFIG_FILE="/home2/adefalco/LiquidBiopsy_Data/WGS/TEST_PIPELINE_VARIANT/Config/ENV_USER_OPTIONS.yaml"
SAMPLE_IN_PARALLEL="5"

#### ALIGNMENT #####
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s 1_pipelineAlign_noPlot.Snakefile \
	--jobs $SAMPLE_IN_PARALLEL --use-conda

#### FRAGMENTOMICS ####
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s 1_pipelineAlign_noPlot.Snakefile \
	-R all_Fragmentomics \
	--jobs $SAMPLE_IN_PARALLEL --use-conda

#### VARIANTS ####
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s 1_pipelineAlign_noPlot.Snakefile \
	-R --until Do_Variant_Calling \
	--jobs $SAMPLE_IN_PARALLEL --use-conda

#### MULTIQC ####
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s 1_pipelineAlign_noPlot.Snakefile \
	-R --until multiqc_all \
	--jobs $SAMPLE_IN_PARALLEL --use-conda
