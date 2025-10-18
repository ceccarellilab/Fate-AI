
setup_environment <- function(config_path = "Config/config.yaml") {
  # Load libraries
  required_pkgs <- c("yaml", "caret", "TCGAbiolinks", "SummarizedExperiment","parallel", 
                     "doParallel", "pgxRpi", "dplyr", "GenomicRanges", "philentropy", "deconvR")
  invisible(lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }))
  
  # Load config file
  config <- read_yaml(config_path)
  list2env(config, envir = .GlobalEnv)
  message("Config file loaded")
  
  # Load Fate-AI functions
  lapply(as.list(list.files(paste0(config$PATH_INITIAL, "R/"), pattern = ".R")), function(x) source(paste0(config$PATH_INITIAL, "R/",x)))
  message("Fate-AI functions loaded")
  
  ### Identify DMRs from TCGA and Methylation Atlas [Fate-AI(+Meth)]
  saveDMRs_fromTCGA(
    CancerTypes = as.character(CLASS_TO_TCGA), 
  )
  message("DMRs from TCGA obtained")

  ### Generate BED Files for Top DMRs [Fate-AI(+Meth)]
  saveBED_TopDMRs(
    ClassTypes = c("Plasma", names(CLASS_TO_TCGA))
  )
  message("BED files of DMRs Generated")
  
  message("✅ Environment successfully initialized.")
}



# # Define mapping of class to TCGA
# CLASS_TO_TCGA <- list(
#   "Colon"       = "TCGA_COAD",
#   "Lung_LUAD"   = "TCGA-LUAD",
#   "Lung_LUSC"   = "TCGA-LUSC",
#   "Breast"      = "TCGA-BRCA",
#   "Prostate"    = "TCGA-PRAD",
#   "Melanoma"    = "TCGA-SKCM",
#   "Mesothelioma" = "TCGA-MESO",
#   "Urothelial"  = "TCGA-BLCA"
# )
# 
# # Define mapping of class to frequency and file path
# CLASS_PARAMS_WGS <- list(
#   Colon          = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C2955.tsv")),
#   Lung           = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4878.tsv")),
#   Prostate       = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4863.tsv")),
#   Urothelial     = list(freq = 22.5,file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4030.tsv")),
#   EwS            = list(freq = 10,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4817.tsv")),
#   Mesothelioma   = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4456.tsv")),
#   Melanoma       = list(freq = 20,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C3224.tsv")),
#   Breast         = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C4872.tsv")),
#   MM             = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C3242.tsv")),
#   Pancreatic     = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/NCIT_C8294.tsv"))
# )

# FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
# PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"
# NUM_THREADS <- 10
# PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
# 
# 
# ALL_BAM_MEDIP_DIR = "/home3/adefalco/Fate-AI/FIGURE/FIGURES_paper/test_script/TEST_BAM_MEDIP"
# SUFFIX_BAM_MEDIP <- ".sorted.bam"
# ALL_BAM_WGS_DIR = "/home2/adefalco/Fate-AI/WGS_alignment/output_folder/BAM"
# SUFFIX_BAM_WGS = "_recal.bam"
# 
# 
# getPathBam <- function(MEDIP = FALSE){
#   ifelse(MEDIP, ALL_BAM_DIR <- ALL_BAM_MEDIP_DIR, ALL_BAM_DIR <- ALL_BAM_WGS_DIR)
#   
#   ALL_BAM_PATH <- list.files(ALL_BAM_DIR, pattern = ".bam.bai", full.names = TRUE)
#   ALL_BAM_PATH <- gsub(".bam.bai", ".bam", ALL_BAM_PATH)
#   ALL_BAM_PATH
# }
# 
# getSamples <- function(MEDIP = FALSE){
#   AllSample <- getPathBam(MEDIP)
#   
#   ifelse(MEDIP, ALL_BAM_DIR <- ALL_BAM_MEDIP_DIR, ALL_BAM_DIR <- ALL_BAM_WGS_DIR)
#   ifelse(MEDIP, SUFFIX_BAM <- SUFFIX_BAM_MEDIP, SUFFIX_BAM <- SUFFIX_BAM_WGS)
#   
#   AllSample <- gsub(ALL_BAM_DIR, "", AllSample)
#   AllSample <- gsub(SUFFIX_BAM, "", AllSample)
#   AllSample <- gsub("/", "", AllSample)
#   AllSample
# }
# 







