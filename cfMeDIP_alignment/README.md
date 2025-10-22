# cfMeDIP Alignment

The entire pipeline is implemented using Snakemake and is executed on the cluster using Slurm.

## 1) Set parameters in Config File (Config/config.yaml): 

### Output & Input Paths

| Parameter | Description | Example |
|-----------|-------------|---------|
| `OUTDIR` | Folder where all mapping outputs (BAM, logs) will be stored | `"output_folder/"` |
| `DATDIR` | Folder containing raw FASTQ sequencing data | `"FASTQ/"` <br>Supports nested structure like `DATDIR/EXP/SAMPTYPE` (e.g., `HOME/DATI/GSE87096/CONTROL`) |
| `PATH_FASTQ` | Path to FASTQ files (can be same as `DATDIR`) | `"FASTQ/"` |

---

### Pipeline Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `NUM_THREADS` | Number of threads to use for parallel processing in each sample | `40` |

---

### Reference Genome & Index

| Parameter | Description | Example |
|-----------|-------------|---------|
| `INDEX` | Path to Bowtie2 genome index | `bcbio/genomes/Hsapiens/hg38/bowtie2` |
| `GENOME_index` | Reference genome name used for alignment | `"hg38"` |

---

### Tool Paths

| Tool | Description | Example |
|------|-------------|---------|
| `SAMTOOLSPATH` | Path to `samtools` binary | `"bcbio/anaconda/bin/samtools"` |
| `BOWTIE2PATH` | Path to `bowtie2` binary | `"bcbio/anaconda/bin/bowtie2"` |
| `FASTQCPATH` | Path to FastQC for quality control | `export PATH="/home/adefalco/cfBreast/FastQC:$PATH"` |
| `CUTADAPTPATH` | Path to `cutadapt` for adapter trimming | `"bcbio/anaconda/bin/cutadapt"` |
| `TRIMGALOREPATH` | Path to `Trim Galore` wrapper | `"tools/TrimGalore-0.6.6/trim_galore"` |

---

### FASTQ Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FASTQ_SUFFIXES` | File name suffixes for paired-end reads | `1_SUFFIX: "R1_001.fastq.gz"` <br> `2_SUFFIX: "R2_001.fastq.gz"` |

---

### Data

| Parameter | Description | Example |
|-----------|-------------|---------|
| `EXP` | List of experimental datasets to process | `batch1` |
| `SAMPTYPE` | List of sample types in each dataset | `batch1: Melanoma, Healthy` |
| `EXCLUDE_SAMPLES` | List of samples to exclude from analysis | *(leave empty if none)* |

---


## 2) Run pipeline: 

```
CONFIG_FILE="Fate-AI/cfMeDIP_alignment/Config/config.yaml"
FIXED_CONDA_PATH="Fate-AI/cfMeDIP_alignment/.snakemake/conda/"
SAMPLE_IN_PARALLEL="4"

snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s cfMedip_split.Snakefile \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
```