# Fragmentomics Analysis for Tumor Evaluation with AI (Fate-AI)

## 1) Install Snakemake and SLURM plugin (tested with snakemake 9.12.0)
```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake snakemake-executor-plugin-slurm
```

## 2) Alignment WGS 
Follow the instructions: [Snakemake pipeline WGS alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/WGS_alignment)

## 3) Alignment cfMeDIP-seq
Follow the instructions: [Snakemake pipeline cfMeDIP alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/cfMeDIP_alignment)


## 4) Fate-AI(+Meth)

### Set Environment Variables

The Fate-AI(+Meth) pipeline requires several R environment variables to define paths, parameters, and genomic references.

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `CLASS_TO_TCGA` | list | Maps local cancer class names to TCGA cohort identifiers. | `Colon = "TCGA_COAD"` |
| `CLASS_PARAMS_WGS` | list | Parameters for each tumor type: frequency of samples (`freq`) and reference file path (`file`). | `Colon = list(freq=25, file=paste0(PATH_INITIAL,"data/progenetix/NCIT_C2955.tsv"))` |
| `FASTA_FILE` | string | Path to the human genome FASTA file (hg38). | `"/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"` |
| `PATH_SAMTOOLS` | string | Path to the `samtools` binary. | `"/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"` |
| `NUM_THREADS` | integer | Number of threads for parallel processing. | `10` |
| `PATH_INITIAL` | string | Root directory of the pipeline project. | `"/home2/adefalco/Fate-AI/"` |

### Fate-AI.1) Compute fragment lengths in each bin (3Mb)
```
AllSample_df <- data.frame(
  Sample = "ICH20", 
  pathBAM_WGS = "raw_data/BAM/WGS/ICH20_recal.bam", 
  pathBAM_MEDIP = "raw_data/BAM/cfMeDIP/IPI01S.sorted.bam", 
  Class = "Urothelial", 
  row.names = "ICH20"
)

saveFragmBIN_fromBam(
        PATH_INITIAL = PATH_INITIAL, 
        sample = AllSample_df$Sample[i], 
        bam = AllSample_df$pathBAM_WGS[i], 
        NUM_THREADS = NUM_THREADS, 
        PATH_SAMTOOLS = PATH_SAMTOOLS, 
        FASTA_FILE = FASTA_FILE, 
        SUFFIX_BAM = gsub(".bam","", SUFFIX_BAM_WGS)
    )
```

### Fate-AI.2) Compute fragment lengths in each bin (3Mb)
```
 saveMetricsBIN(
        PATH_INITIAL = PATH_INITIAL, 
        sample = AllSample_df$Sample[i],
        NUM_THREADS = NUM_THREADS
    )
```

### Meth.1) Identify DMRs from TCGA and Methylation Atlas (+Meth)
```
saveDMRs_fromTCGA(
  PATH_INITIAL = PATH_INITIAL, 
  CancerTypes = as.character(CLASS_TO_TCGA), 
  NUM_THREADS = NUM_THREADS
)
```
### Meth.2) Generate BED Files for Top DMRs (+Meth)

```
saveBED_TopDMRs(
  PATH_INITIAL = PATH_INITIAL, 
  ClassTypes = c("Plasma", names(CLASS_TO_TCGA))
)
```

### Meth.3) Generate BED Files for Top DMRs (+Meth)

```
saveBED_TopDMRs(
  PATH_INITIAL = PATH_INITIAL, 
  ClassTypes = c("Plasma", names(CLASS_TO_TCGA))
)
```

### Meth.4) Get Coverage on DMRs for Each Sample (+Meth)

```
saveCoverageDMRs_fromBam(
    PATH_INITIAL = PATH_INITIAL, 
    sample = AllSample_df$Sample[i],
    bam = AllSample_df$pathBAM_MEDIP[i],
    FASTA_FILE = FASTA_FILE,
    PATH_SAMTOOLS = PATH_SAMTOOLS,
    ClassTypes = c(
      AllSample_df$Class[i]
    )
  )
```


### Fate-AI(+Meth).1) Get features
```
feat_WGS <- getFeatureBasedOnCNV(
    AllSample_df$Sample, 
    PATH_INITIAL = PATH_INITIAL, 
    CLASS_CNV = AllSample_df$Class[1], 
    NUM_THREADS = NUM_THREADS
)
feat_cfmedip <- getFeature_cfMeDIP(
    AllSample_df$Sample,
    PATH_INITIAL = PATH_INITIAL,
    CLASS = AllSample_df$Class[1]
)
feat_mtx <- cbind(feat_WGS, feat_cfmedip)
```

### Fate-AI(+Meth).2) Get prediction

```
TEST_INDEX <- NULL

if(is.null(TEST_INDEX)){ #"Matched-cohort"
    prediction <- classifyMATRIX(
        feat_mtx, 
        classes = AllSample_SUB$Class, 
        class1 = "Healthy", 
        class2 = CLASS, 
        method = METHOD_CLASSIFIER
    )
} else { #"Cross-cohort"
    prediction <- classifyMATRIX(
        feat_mtx, 
        classes = AllSample_SUB$Class, 
        class1 = "Healthy", 
        class2 = CLASS, 
        testIND = TEST_INDEX, 
        method = METHOD_CLASSIFIER
    )
}
```