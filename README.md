<div>
<img src="Fate-AI.png" align="right" width="200">
</div>

# Fragmentomics Analysis for Tumor Evaluation with AI (Fate-AI)

Preprint: [Knowledge-informed multimodal cfDNA analysis improves sensitivity and generalization in cancer detection](https://www.biorxiv.org/content/10.1101/2025.10.20.683167v1)

## 1) Install Fate-AI and Snakemake

Install snakemake and Slurm executor plugin for pipeline alignment (tested with snakemake 9.12.0), and clone the Fate-AI repository:
```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake snakemake-executor-plugin-slurm
git clone https://github.com/ceccarellilab/Fate-AI.git
```

Install Fate-AI in R:
```
library(devtools)
install_github("ceccarellilab/Fate-AI")
```

## 2) Align lpWGS samples and perform GC correction
To align lpWGS samples or perform GC correction on aligned files, follow the instructions in the [Snakemake pipeline WGS alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/WGS_alignment)

## 3) Alignment cfMeDIP-seq
To align cfMeDIP-seq samples, follow the instructions: [Snakemake pipeline cfMeDIP alignment](https://github.com/ceccarellilab/Fate-AI/tree/main/cfMeDIP_alignment)

## 4) Fate-AI(+Meth)

### 4.1) Set parameters in Config File (Config/config.yaml): 

### Paths and General Settings

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `PATH_INITIAL` | string | Root directory of the pipeline project. | `"Fate-AI/"` |
| `FASTA_FILE` | string | Path to the human genome FASTA file (hg38). | `"bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"` |
| `PATH_SAMTOOLS` | string | Path to the `samtools` binary. | `"bin/samtools"` |
| `BED_TOOLS_DIR` | Path to `bedtools` binary | `"bcbio/anaconda/bin/bedtools"` |
| `NUM_THREADS` | integer | Number of threads for parallel processing. | `10` |

### Primary Tumor Mapping

| Variable | Type | Description | Example |
|----------|------|-------------|---------|
| `CLASS_PARAMS_WGS` | list | Parameters for each cancer type: CNV frequency (`freq`) and NCIT code (`NCIT`). [Fate-AI]  | `Colon: freq: 25 NCIT: "C2955"` |
| `CLASS_TO_TCGA` | list | Maps cancer type to TCGA identifiers. [Fate-AI(+Meth)] | `Urothelial : "TCGA_BLCA"` |

### 4.2) Setup environment 
```
library(FateAI)
setup_environment()
```
### 4.3) Compute fragment lengths and metrics in each bin (3Mb) [Fate-AI]
```
AllSample_df <- data.frame(
  Sample = "ICH20", 
  pathBAM_WGS = "raw_data/BAM/WGS/ICH20_recal.bam", 
  pathBAM_MEDIP = "raw_data/BAM/cfMeDIP/IPI01S.sorted.bam", 
  Class = "Urothelial", 
  row.names = "ICH20"
)

saveMetricsBIN(
        sample = AllSample_df$Sample[i], 
        bam = AllSample_df$pathBAM_WGS[i]
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
        class2 = "Cancer", 
        method = METHOD_CLASSIFIER
    )
} else { #"Cross-cohort"
    prediction <- classifyMATRIX(
        feat_mtx, 
        classes = AllSample_SUB$Class, 
        class1 = "Healthy", 
        class2 = "Cancer", 
        testIND = TEST_INDEX, 
        method = METHOD_CLASSIFIER
    )
}
```



## Citation

> @article {De Falco2025.10.20.683167,  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; author = {De Falco, Antonio and Grisolia, Piera and Giuffrida, Raffaella and Iannarone, Clara and Graziano, Cinzia and Scrima, Marianna and Cossu, Alessia Maria and Tufano, Rossella and Yow, Maria Vanessa and Brown, Chloe Marissa and Bajaj, Palak and Bocchetti, Marco and Nuzzo, Pier Vitale and Morgillo, Floriana and Caraglia, Francesco and Della Corte, Carminia Maria and Di Guida, Gaetano and Troiani, Teresa and Ciardiello, Fortunato and Rizzo, Maria Rosaria and Fiorelli, Alfonso and Giordano, Noemi Maria and Arcaniolo, Davide and Della Rosa, Giampiero and Desio, Marco and Covre, Alessia and Di Giacomo, Annamaria and Calabrò, Luana and Fontana, Paolo and Mare, Mariza and Landgren, Ola and Green, Damian and Lesokhin, Alex and Maio, Michele and Coffey, David and Merchant, Nipun and Datta, Jashodeep and Forte, Stefano and Caraglia, Michele and Ceccarelli, Michele},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; title = {Knowledge-informed multimodal cfDNA analysis improves sensitivity and generalization in cancer detection},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; year = {2025},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; month = {10},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; elocation-id = {2025.10.20.683167},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; doi = {10.1101/2025.10.20.683167},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; URL = {https://www.biorxiv.org/content/early/2025/10/21/2025.10.20.683167},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; eprint = {https://www.biorxiv.org/content/early/2025/10/21/2025.10.20.683167.full.pdf},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; publisher = {Cold Spring Harbor Laboratory},  
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; journal = {bioRxiv}  
> }
