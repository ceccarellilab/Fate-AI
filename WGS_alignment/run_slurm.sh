CONFIG_FILE="/home2/adefalco/LiquidBiopsy_Data/WGS/TEST_PIPELINE_VARIANT/Config/ENV_USER_OPTIONS.yaml"
FIXED_CONDA_PATH="/home3/adefalco/Fate-AI/output_folder/.snakemake/conda/"
SAMPLE_IN_PARALLEL="4"

#### ALIGNMENT #####
#snakemake --configfile $CONFIG_FILE \
#	--executor slurm --workflow-profile profiles/default/ \
#	-s 1_pipelineAlign_noPlot.Snakefile \
#	--jobs $SAMPLE_IN_PARALLEL \
#	--use-conda --conda-prefix $FIXED_CONDA_PATH

#### FRAGMENTOMICS ####
/scratch/adefalco/miniconda3/envs/snakemake_last/bin/snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s 1_pipelineAlign_noPlot.Snakefile \
	-R all_Fragmentomics \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH

#### VARIANTS ####
#snakemake --configfile $CONFIG_FILE \
#	--executor slurm --workflow-profile profiles/default/ \
#	-s 1_pipelineAlign_noPlot.Snakefile \
#	-R --until Do_Variant_Calling \
#	--jobs $SAMPLE_IN_PARALLEL \
#	--use-conda --conda-prefix $FIXED_CONDA_PATH

#### MULTIQC ####
#snakemake --configfile $CONFIG_FILE \
#	--executor slurm --workflow-profile profiles/default/ \
#	-s 1_pipelineAlign_noPlot.Snakefile \
#	-R --until multiqc_all \
#	--jobs $SAMPLE_IN_PARALLEL \
#	--use-conda --conda-prefix $FIXED_CONDA_PATH
