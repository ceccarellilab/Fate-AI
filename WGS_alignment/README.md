# lpWGS alignment

The entire pipeline is implemented using Snakemake and is executed on the cluster using Slurm.

The pipeline can start from either raw sequencing FASTQ files or pre-aligned BAM files.
- Starting from ***FASTQ***: allows end-to-end processing, including alignment, quality filtering, deduplication and recalibration using Sentieon, and finally GC correction.
- Starting from ***BAM***: files skips the alignment step and proceeds directly to GC correction.

## 1) Set parameters in Config File (Config/ENV_USER_OPTIONS.yaml): 

### Output & Input Paths

| Parameter | Description | Example |
|-----------|-------------|---------|
| `WORKPATH` | Folder where all output files will be stored | `"output_folder/"` |
| `PATH_FASTQ` | Folder containing raw FASTQ files (if starting from sequencing reads) | `"FASTQ/"` |
| `PATH_BAM` | Folder containing pre-aligned BAM files (leave empty if starting from FASTQ) | `"BAM_FILE"` |

---

### FASTQ Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `PAIRED` | Indicates if sequencing data is paired-end (`True`) or single-end (`False`) | `True` |
| `FASTQ_SUFFIXES` | File name suffixes for read 1 and read 2 | `1_SUFFIX: "R1_001.fastq.gz"`<br>`2_SUFFIX: "R2_001.fastq.gz"` |
| `PLATFORM` | Sequencing platform | `"ILLUMINA"` |

---

### BAM Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `BAM_SUFFIX` | Suffix for BAM files | `"_recal.bam"` |
| `NUM_THREADS` | Number of threads to use for parallel processing in each sample | `40` |
| `MAX_FRAG_LENGHT` | Maximum fragment length to consider for fragment-based analysis | `550` |
| `MAPQ` | Minimum mapping quality threshold for reads | `30` |
| `BIN_SIZE` | Bin size for genome-wide metrics (e.g., copy number or fragmentation) | `3000000` |

---

### Reference Data Files

| Parameter | Description | Example |
|-----------|-------------|---------|
| `FASTA` | Reference genome FASTA file | `"bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"` |
| `DBSNP` | Known SNP VCF file for variant calling and annotation | `"bcbio/genomes/Hsapiens/hg38/variation/dbsnp-150.vcf.gz"` |
| `KNOWN_SITES` | Known variant sites for base recalibration | `MILLS_INDELS: "…Mills_and_1000G_gold_standard.indels.vcf.gz"`<br>`1000G_INDELS: "…1000G_phase1.snps.high_confidence.vcf.gz"` |

---

### Sentieon Configuration 

| Parameter | Description | Example |
|-----------|-------------|---------|
| `SENTIEON_MODULE` | Name of the cluster module for Sentieon (leave empty if not using modules) | `"sentieon_202503"` |
| `SENTIEON_PATH` | Path to Sentieon binaries if modules are not used | `"/storage/qnap_vol1/SHARED/NGSTOOLS/SENTIEON/sentieon-genomics-202503"` |
| `SENTIEON_LICENSE` | License server for Sentieon WGS tools | `"XXX.XXX.X.XXX:XXXX"` |
| `ML_MODEL_N` | Path to Sentieon machine learning model bundle (Optional) | `"SentieonIlluminaWGS2.2.bundle"` |
| `TNSCOPE_FILTER` | Path to TNScope filtering script (Optional) | `"tnscope_filter.py"` |

---

### Tool Paths

| Tool | Description | Example |
|------|-------------|---------|
| `BCF_DIR` | Path to `bcftools` binary | `"bcbio/anaconda/bin/bcftools"` |
| `BGZIP_DIR` | Path to `bgzip` binary | `"bcbio/anaconda/bin/bgzip"` |
| `SNPEFF_DIR` | Path to `snpEff` binary | `"bcbio/anaconda/bin/snpEff"` |
| `SAMTOOLS_DIR` | Path to `samtools` binary | `"bcbio/anaconda/bin/samtools"` |


---


## 2) Run pipeline: 

```
CONFIG_FILE="Fate-AI/WGS_alignment/Config/ENV_USER_OPTIONS.yaml"
FIXED_CONDA_PATH="Fate-AI/WGS_alignment/output_folder/.snakemake/conda/"
SNAKEFILE="Fate-AI/WGS_alignment/pipelineAlignWGS.Snakefile"
SAMPLE_IN_PARALLEL="4"

conda activate snakemake
snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s $SNAKEFILE \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
```

## 3) Run variant calling (optional):  
```
snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s $SNAKEFILE \
	-R --until Do_Variant_Calling \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
```

## 4) Run multiqc (optional):  
```
snakemake --configfile $CONFIG_FILE \
	--executor slurm --workflow-profile profiles/default/ \
	-s $SNAKEFILE \
	-R --until multiqc_all \
	--jobs $SAMPLE_IN_PARALLEL \
	--use-conda --conda-prefix $FIXED_CONDA_PATH
```
