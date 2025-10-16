##### WORKFLOW cfMEDIP#####

PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))

#Get DMRs from TCGA and Methylation Atlas Deconvolution
saveDMRs_fromTCGA(PATH_INITIAL = PATH_INITIAL, CancerTypes = as.character(CLASS_TO_TCGA), NUM_THREADS = NUM_THREADS)

#Generate BED files of DMRs
saveBED_TopDMRs(PATH_INITIAL = PATH_INITIAL, ClassTypes = c("Plasma", names(CLASS_TO_TCGA)))

AllSample_df <- data.frame(sample = getSamples(MEDIP = T), pathBAM = getPathBam(MEDIP = T), row.names = getSamples(MEDIP = T))

(load("~/home3/Fate-AI/FIGURE/FIGURES_paper/data/AllSample_ALL_with_coverage_MEDIP.RData"))
AllSample_ALL <- AllSample_ALL["IPH20",]

AllSample_df <- data.frame(sample = AllSample_ALL$IC_code, pathBAM_WGS = AllSample_ALL$PATH_BAM_WGS, pathBAM_MEDIP = AllSample_ALL$PATH_BAM_MEDIP, Class = AllSample_ALL$Class, row.names = AllSample_ALL$IC_code)


lapply(1:nrow(AllSample_df), function(i){
  
  #Get coverage on DMRs for each sample
  saveCoverageDMRs_fromBam(PATH_INITIAL = PATH_INITIAL, 
                           sample = AllSample_df$sample[i],
                           bam = AllSample_df$pathBAM_MEDIP[i],
                           FASTA_FILE = FASTA_FILE,
                           PATH_SAMTOOLS = PATH_SAMTOOLS,
                           #AllSample_df$Class)
                           ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma", "Lung_SHARED"))
})

#Get features (cfMeDIP)
feat_cfmedip <- getFeature_cfMeDIP(AllSample_df$sample,
                                   PATH_INITIAL = PATH_INITIAL,
                                   #CLASS = AllSample_df$Class)
                                   CLASS = c("Colon", "Prostate", "Breast", "Lung", "Mesotelioma", "Melanoma", "Urothelial"))


##### WORKFLOW WGS #####

PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
lapply(as.list(list.files(paste0(PATH_INITIAL, "scripts/"), pattern = ".R")), function(x) source(paste0(PATH_INITIAL, "scripts/",x)))

AllSample_df <- data.frame(sample = getSamples(), pathBAM = getPathBam(), row.names = getSamples())

lapply(1:nrow(AllSample_df), function(i){

saveFragmBIN_fromBam(PATH_INITIAL = PATH_INITIAL, 
                     sample = AllSample_df$sample[i], 
                     bam = AllSample_df$pathBAM[i], 
                     NUM_THREADS = NUM_THREADS, 
                     PATH_SAMTOOLS = PATH_SAMTOOLS, 
                     FASTA_FILE = FASTA_FILE, 
                     SUFFIX_BAM = gsub(".bam","", SUFFIX_BAM_WGS))
  
saveMetricsBIN(PATH_INITIAL = PATH_INITIAL, 
                 sample = AllSample_df$sample[i],
                 NUM_THREADS = NUM_THREADS)

})

feat_WGS <- getFeatureBasedOnCNV(AllSample_df$sample, PATH_INITIAL = PATH_INITIAL, 
                                 CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                 NUM_THREADS = NUM_THREADS)


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
  medip_mtx <- getFeatMEDIP_last_COUNT(AllSample, )
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
  

    


