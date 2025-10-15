##### WORKFLOW cfMEDIP#####

PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))

saveDMRs_fromTCGA(PATH_INITIAL = PATH_INITIAL, CancerTypes = as.character(CLASS_TO_TCGA), NUM_THREADS = NUM_THREADS)

saveBED_TopDMRs(PATH_INITIAL = PATH_INITIAL, ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma"))

saveCoverageDMRs_fromBam(PATH_INITIAL = PATH_INITIAL, 
                         ALL_BAM_MEDIP_PATH = ALL_BAM_MEDIP_PATH, 
                         FASTA_FILE = FASTA_FILE,
                         PATH_SAMTOOLS = PATH_SAMTOOLS,
                         ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma", "Lung_SHARED"),
                         SUFFIX_BAM = SUFFIX_BAM, 
                         NUM_THREADS = NUM_THREADS)


merge_Samples_cfMEDIP_Counts(PATH_INITIAL = PATH_INITIAL)

feat_cfmedip <- getFeature_cfMeDIP(PATH_INITIAL = PATH_INITIAL, ALL_BAM_MEDIP_PATH = ALL_BAM_MEDIP_PATH, SUFFIX_BAM = SUFFIX_BAM,
                                   ClassTypes = c("Plasma", "Colon", "Prostate", "Breast", "Lung", "Mesotelioma", "Melanoma", "Urothelial"))


##### WORKFLOW WGS #####

PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))


# remotes::install_github("progenetix/pgxRpi")
# library("pgxRpi")
# frequency <- pgxLoader(type="cnv_frequency", output ='pgxfreq',
#                        filters=c("NCIT:C3224"))
# pgxFreqplot(frequency)

AllSample_df <- data.frame(sample = getSamples(), pathBAM = getPathBam(), row.names = getSamples())

EXAMPLE_SAMPLE <- AllSample_df$sample[1]
SAMPLE <- AllSample_df[EXAMPLE_SAMPLE,]$sample
BAM <- AllSample_df[EXAMPLE_SAMPLE,]$pathBAM
saveFragmBIN_fromBam(PATH_INITIAL = PATH_INITIAL, sample = SAMPLE, bam = BAM, NUM_THREADS = NUM_THREADS, PATH_SAMTOOLS = PATH_SAMTOOLS, FASTA_FILE = FASTA_FILE, SUFFIX_BAM = gsub(".bam","", SUFFIX_BAM))

saveMetricsBIN(PATH_INITIAL = PATH_INITIAL, 
               sample = SAMPLE,
               NUM_THREADS = NUM_THREADS)



feat_WGS <- getFeatureBasedOnCNV(AllSample, PATH_INITIAL = PATH_INITIAL, 
                                 CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                 NUM_THREADS = 30)


#### Feature WGS and Medip ####

lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))

CLASS <- "Colon"
METHOD_CLASSIFIER <- "glmnet"
MODEL <- "Fate-AI(+Meth)"
NUM_THREADS <- 10

if(MODEL == "Fate-AI"){
  
  feat_mtx <- getFeatureBasedOnCNV(AllSample, PATH_INITIAL = PATH_INITIAL, 
                                   CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                   NUM_THREADS = NUM_THREADS)
  METHOD_CLASSIFIER <- "glmnet"
  
}else if(MODEL == "Fate-AI(+Meth)"){
  
  feat_WGS <- getFeatureBasedOnCNV(AllSample, PATH_INITIAL = PATH_INITIAL, 
                                   CLASS_CNV = CLASS, 
                                   NUM_THREADS = NUM_THREADS)
  medip_mtx <- getFeatMEDIP_last_COUNT(AllSample)
  feat_mtx <- cbind(feat_mtx, medip_mtx)

}

METHOD_CLASSIFIER <- "glmnet"

feat_mtx <- normalizeMatrixDataset(feat_mtx, AllSample)

file_name = paste(MODEL, CLASS, METHOD_CLASSIFIER, sep = "_")

file_name <- paste0(file_name, ".RData")

if(is.null(TEST_INDEX)){
  prediction <- classifyMATRIX(feat_mtx, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, method = METHOD_CLASSIFIER)
  TYPE = "Matched-cohort"
}else{
  prediction <- classifyMATRIX(feat_mtx, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, testIND = TEST_INDEX, method = METHOD_CLASSIFIER)
  TYPE = "Cross-cohort"
}

AUC <- getAUC(prediction, label0 = CLASS)

resAUCs <- data.frame(AUC = AUC, TYPE = TYPE, METHOD_CLASSIFIER = METHOD_CLASSIFIER, MODEL = MODEL)

list(prediction, resAUCs)
  

    


