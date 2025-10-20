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

### 4.1) Set parameters in Config File (Config/config.yaml): 

### Paths and General Settings

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `PATH_INITIAL` | string | Root directory of the pipeline project. | `"Fate-AI/"` |
| `FASTA_FILE` | string | Path to the human genome FASTA file (hg38). | `"bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"` |
| `PATH_SAMTOOLS` | string | Path to the `samtools` binary. | `"bin/samtools"` |
| `NUM_THREADS` | integer | Number of threads for parallel processing. | `10` |

### Primary Tumor Mapping

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `CLASS_PARAMS_WGS` | list | Parameters for each cancer type: CNV frequency (`freq`) and NCIT code (`NCIT`). [Fate-AI]  | `Colon: freq: 25 NCIT: "C2955"` |
| `CLASS_TO_TCGA` | list | Maps cancer type to TCGA identifiers. [Fate-AI(+Meth)] | `Urothelial : "TCGA_BLCA"` |

### 4.2) Compute fragment lengths in each bin (3Mb) [Fate-AI]
```
AllSample_df <- data.frame(
  Sample = "ICH20", 
  pathBAM_WGS = "raw_data/BAM/WGS/ICH20_recal.bam", 
  pathBAM_MEDIP = "raw_data/BAM/cfMeDIP/IPI01S.sorted.bam", 
  Class = "Urothelial", 
  row.names = "ICH20"
)

saveFragmBIN_fromBam(
        sample = AllSample_df$Sample[i], 
        bam = AllSample_df$pathBAM_WGS[i]
    )
```

### 4.3) Compute metrics in each bin (3Mb) [Fate-AI]
```
 saveMetricsBIN(
        sample = AllSample_df$Sample[i],
    )
```

### 4.4) Get Coverage on DMRs for Each Sample [Fate-AI(+Meth)]

```
saveCoverageDMRs_fromBam(
    sample = AllSample_df$Sample[i],
    bam = AllSample_df$pathBAM_MEDIP[i],
    ClassTypes = c(
      AllSample_df$Class[i]
    )
  )
```


###  4.5) Get features lpWGS [Fate-AI]
```
feat_WGS <- getFeatureBasedOnCNV(
    AllSample_df$Sample, 
    CLASS_CNV = AllSample_df$Class[1]
)
```

###  4.6) Get features cfMeDIP [Fate-AI(+Meth)]
```
feat_cfmedip <- getFeature_cfMeDIP(
    AllSample_df$Sample,
    CLASS = AllSample_df$Class[1]
)

```

### 4.7) Get prediction [Fate-AI(+Meth)]

```
feat_mtx <- cbind(feat_WGS, feat_cfmedip)
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
