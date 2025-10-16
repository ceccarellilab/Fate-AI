# Define mapping of class to TCGA
CLASS_TO_TCGA <- list(
  "Colon"       = "TCGA_COAD",
  "Lung_LUAD"   = "TCGA-LUAD",
  "Lung_LUSC"   = "TCGA-LUSC",
  "Breast"      = "TCGA-BRCA",
  "Prostate"    = "TCGA-PRAD",
  "Melanoma"    = "TCGA-SKCM",
  "Mesothelioma" = "TCGA-MESO",
  "Urothelial"  = "TCGA-BLCA"
)

# remotes::install_github("progenetix/pgxRpi")
# library("pgxRpi")
# frequency <- pgxLoader(type="cnv_frequency", output ='pgxfreq',
#                        filters=c("NCIT:C3224"))
# pgxFreqplot(frequency)

# Define mapping of class to frequency and file path
CLASS_PARAMS_WGS <- list(
  Colon          = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Colorectal_Carcinoma_NCIT_C2955.tsv")),
  Lung           = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Lung_Carcinoma_NCIT_C4878.tsv")),
  Prostate       = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Prostate_Carcinoma_NCIT_C4863.tsv")),
  Urothelial     = list(freq = 22.5,file = paste0(PATH_INITIAL, "data/progenetix/Urothelial_Carcinoma_NCIT_C4030.tsv")),
  EwS            = list(freq = 10,  file = paste0(PATH_INITIAL, "data/progenetix/Ewing_Sarcoma_NCIT_C4817.tsv")),
  Mesothelioma   = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Malignant_Mesothelioma_NCIT_C4456.tsv")),
  Melanoma       = list(freq = 20,  file = paste0(PATH_INITIAL, "data/progenetix/Melanoma_NCIT_C3224.tsv")),
  Breast         = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Breast_Carcinoma_NCIT_C4872.tsv")),
  MM             = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Multiple_Myeloma_NCIT_C3242.tsv")),
  Pancreatic     = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Pancreatic_Adenocarcinoma_NCIT_C8294.tsv"))
)

FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"
NUM_THREADS <- 10
PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
ALL_BAM_MEDIP_DIR = "/home3/adefalco/Fate-AI/FIGURE/FIGURES_paper/test_script/TEST_BAM_MEDIP"
SUFFIX_BAM_MEDIP <- ".sorted.bam"

ALL_BAM_WGS_DIR = "/home2/adefalco/Fate-AI/WGS_alignment/output_folder/BAM"
SUFFIX_BAM_WGS = "_recal.bam"

getPathBam <- function(MEDIP = FALSE){
  ifelse(MEDIP, ALL_BAM_DIR <- ALL_BAM_MEDIP_DIR, ALL_BAM_DIR <- ALL_BAM_WGS_DIR)
  
  ALL_BAM_PATH <- list.files(ALL_BAM_DIR, pattern = ".bam.bai", full.names = TRUE)
  ALL_BAM_PATH <- gsub(".bam.bai", ".bam", ALL_BAM_PATH)
  ALL_BAM_PATH
}

getSamples <- function(MEDIP = FALSE){
  AllSample <- getPathBam(MEDIP)
  
  ifelse(MEDIP, ALL_BAM_DIR <- ALL_BAM_MEDIP_DIR, ALL_BAM_DIR <- ALL_BAM_WGS_DIR)
  ifelse(MEDIP, SUFFIX_BAM <- SUFFIX_BAM_MEDIP, SUFFIX_BAM <- SUFFIX_BAM_WGS)
  
  AllSample <- gsub(ALL_BAM_DIR, "", AllSample)
  AllSample <- gsub(SUFFIX_BAM, "", AllSample)
  AllSample <- gsub("/", "", AllSample)
  AllSample
}





