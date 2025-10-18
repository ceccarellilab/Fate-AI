############################################# cfMeDIP-seq ##########################################################################

# 1) Prepare the necessary cfMeDIP-seq files, identify DMRs, and generate BED files.

#PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
#lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))

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
  
  # Load Fate-AI functions
  lapply(as.list(list.files(paste0(config$PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(config$PATH_INITIAL, "scripts/",x)))
  
  message("✅ Environment successfully initialized.")
}

setwd("/home2/adefalco/Fate-AI/")
setup_environment()

#Get DMRs from TCGA and Methylation Atlas Deconvolution
saveDMRs_fromTCGA(CancerTypes = as.character(CLASS_TO_TCGA), NUM_THREADS = NUM_THREADS)

#Generate BED files of DMRs
saveBED_TopDMRs(ClassTypes = c("Plasma", names(CLASS_TO_TCGA)))

# 2) Get coverage on DMRs for each sample cfMeDIP-seq

#Data.Frame samples
AllSample_df <- data.frame(Sample = "ICH20", 
                           pathBAM_WGS = "/home/adefalco/ctDNAanalysis/PipelineAling/SnakePipeline/AllData/ICH20_recal.bam", 
                           pathBAM_MEDIP = "/home3/adefalco/Fate-AI/ScriptMedip/BAM_MEDIP/IPI01S.sorted.bam", 
                           Class = "Urothelial", 
                           row.names = "ICH20")

lapply(1:nrow(AllSample_df), function(i){
  
  #Get coverage on DMRs
  saveCoverageDMRs_fromBam(PATH_INITIAL = PATH_INITIAL, 
                           sample = AllSample_df$Sample[i],
                           bam = AllSample_df$pathBAM_MEDIP[i],
                           FASTA_FILE = FASTA_FILE,
                           PATH_SAMTOOLS = PATH_SAMTOOLS,
                           ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma", "Lung_SHARED"))
})

############################################# lpWGS ##########################################################################

lapply(1:nrow(AllSample_df), function(i){

# Get fragment lenght in each bin (3MB)
saveFragmBIN_fromBam(PATH_INITIAL = PATH_INITIAL, 
                     sample = AllSample_df$Sample[i], 
                     bam = AllSample_df$pathBAM_WGS[i], 
                     NUM_THREADS = NUM_THREADS, 
                     PATH_SAMTOOLS = PATH_SAMTOOLS, 
                     FASTA_FILE = FASTA_FILE, 
                     SUFFIX_BAM = gsub(".bam","", SUFFIX_BAM_WGS))

# Get metrics each bin (3MB)  
saveMetricsBIN(PATH_INITIAL = PATH_INITIAL, 
                 sample = AllSample_df$Sample[i],
                 NUM_THREADS = NUM_THREADS)

})

############################################# Fate-AI(+Meth) ##########################################################################

METHOD_CLASSIFIER <- "glmnet"
MODEL <- "Fate-AI(+Meth)"
NUM_THREADS <- 10

CLASS_CNV_METH <- "Urothelial"

if(MODEL == "Fate-AI"){
  
  #Features WGS
  feat_mtx <- getFeatureBasedOnCNV(AllSample$Sample, PATH_INITIAL = PATH_INITIAL, 
                                   CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                   NUM_THREADS = NUM_THREADS)
  METHOD_CLASSIFIER <- "glmnet"
  
}else if(MODEL == "Fate-AI(+Meth)"){
  
  #Features WGS
  feat_WGS <- getFeatureBasedOnCNV(AllSample_df$Sample, PATH_INITIAL = PATH_INITIAL, 
                                   CLASS_CNV = AllSample_df$Class[1], 
                                   NUM_THREADS = NUM_THREADS)
  #Features cfMeDIP-seq
  feat_cfmedip <- getFeature_cfMeDIP(AllSample_df$Sample,
                                     PATH_INITIAL = PATH_INITIAL,
                                     CLASS = AllSample_df$Class)
  
  feat_mtx <- cbind(feat_WGS, feat_cfmedip)

}

METHOD_CLASSIFIER <- "glmnet"

feat_mtx <- normalizeMatrix(feat_mtx)

TEST_INDEX <- NULL

if(is.null(TEST_INDEX)){
  prediction <- classifyMATRIX(feat_mtx, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, method = METHOD_CLASSIFIER)
  TYPE = "Matched-cohort"
}else{
  prediction <- classifyMATRIX(feat_mtx, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, testIND = TEST_INDEX, method = METHOD_CLASSIFIER)
  TYPE = "Cross-cohort"
}

prediction
  

    


