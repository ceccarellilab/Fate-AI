# cfMeDIP Alignment

```
CONFIG_FILE="/home2/adefalco/LiquidBiopsy_Data/cfMedip/pipeline/config_Mel.yaml"
FIXED_CONDA_PATH="/home2/adefalco/Fate-AI/output_folder/.snakemake/conda/"
SAMPLE_IN_PARALLEL="4"
```

## ALIGNMENT
```
snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s cfMedip_split.Snakefile \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
```