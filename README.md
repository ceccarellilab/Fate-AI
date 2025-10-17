# Fragmentomics Analysis for Tumor Evaluation with AI (Fate-AI)

## 1) Install Snakemake and SLURM plugin (tested with snakemake 9.12.0)
```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake snakemake-executor-plugin-slurm
conda activate snakemake
```

## 2) Alignment WGS 
Follow the instructions: [Snakemake pipeline WGS alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/WGS_alignment)

## 3) Alignment cfMeDIP-seq
Follow the instructions: [Snakemake pipeline cfMeDIP alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/cfMeDIP_alignment)

