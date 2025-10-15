##### WORKFLOW cfMEDIP#####

FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"
NUM_THREADS <- 10
INITIAL_PATH <- "/home2/adefalco/Fate-AI"
ALL_BAM_MEDIP_PATH = "/home3/adefalco/Fate-AI/FIGURE/FIGURES_paper/test_script/TEST_BAM_MEDIP"
SUFFIX_BAM = ".sorted.bam"

#TCGA_DIR <- "data/tcga/"
#BED_DIR <- "data/bed_dmrs"
#N_TOP <- 3000

source(paste0(INITIAL_PATH, "scripts/cfMeDIP_features.R"))

CLASS_TO_TCGA <- list(
  "Hyper_Colon"       = "TCGA_COAD",
  "Hyper_Lung_LUAD"   = "TCGA-LUAD",
  "Hyper_Lung_LUSC"   = "TCGA-LUSC",
  "Hyper_Breast"      = "TCGA-BRCA",
  "Hyper_Prostate"    = "TCGA-PRAD",
  "Hyper_Melanoma"    = "TCGA-SKCM",
  "Hyper_Mesotelioma" = "TCGA-MESO",
  "Hyper_Urothelial"  = "TCGA-BLCA"#,
  #"Hyper_Gliobastoma" = "TCGA-GBM"
)


saveDMRs_fromTCGA(INITIAL_PATH = INITIAL_PATH, CancerTypes = as.character(CLASS_TO_TCGA), NUM_THREADS = NUM_THREADS)

saveBED_TopDMRs(INITIAL_PATH = INITIAL_PATH, ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma"))

saveCoverageDMRs_fromBam(INITIAL_PATH = INITIAL_PATH, 
                         ALL_BAM_MEDIP_PATH = ALL_BAM_MEDIP_PATH, 
                         FASTA_FILE = FASTA_FILE,
                         PATH_SAMTOOLS = PATH_SAMTOOLS,
                         ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma", "Lung_SHARED"),
                         SUFFIX_BAM = SUFFIX_BAM, 
                         NUM_THREADS = NUM_THREADS)


merge_Samples_cfMEDIP_Counts(INITIAL_PATH = INITIAL_PATH)

feat_cfmedip <- getFeature_cfMeDIP(INITIAL_PATH = INITIAL_PATH, ALL_BAM_MEDIP_PATH = ALL_BAM_MEDIP_PATH, SUFFIX_BAM = SUFFIX_BAM,
                                   ClassTypes = c("Plasma", "Colon", "Prostate", "Breast", "Lung", "Mesotelioma", "Melanoma", "Urothelial"))


##### WORKFLOW WGS #####
FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"
NUM_THREADS <- 10
PATH_INITIAL <- "/home2/adefalco/Fate-AI/"
ALL_BAM_WGS_DIR = "/home2/adefalco/Fate-AI/WGS_alignment/output_folder/BAM/"
SUFFIX_BAM = "_recal.bam"
source(paste0(PATH_INITIAL, "scripts/lpWGS_features.R"))

# remotes::install_github("progenetix/pgxRpi")
# library("pgxRpi")
# frequency <- pgxLoader(type="cnv_frequency", output ='pgxfreq',
#                        filters=c("NCIT:C3224"))
# pgxFreqplot(frequency)


ALL_BAM_WGS <- list.files(ALL_BAM_WGS_DIR, pattern = ".bam.bai", full.names = TRUE)
ALL_BAM_WGS <- gsub(".bam.bai", ".bam", ALL_BAM_WGS)
AllSample <- ALL_BAM_WGS
AllSample <- gsub(ALL_BAM_WGS_DIR, "", AllSample)
AllSample <- gsub(SUFFIX_BAM, "", AllSample)
AllSample <- gsub("/", "", AllSample)

AllSample_df <- data.frame(sample = AllSample, pathBAM = ALL_BAM_WGS, row.names = AllSample)

EXAMPLE_SAMPLE <- "LB-CRC-32-P-02"
SAMPLE <- AllSample_df[EXAMPLE_SAMPLE,]$sample
BAM <- AllSample_df[EXAMPLE_SAMPLE,]$pathBAM
saveFragmBIN_fromBam(PATH_INITIAL = PATH_INITIAL, sample = SAMPLE, bam = BAM, NUM_THREADS = NUM_THREADS, PATH_SAMTOOLS = PATH_SAMTOOLS, FASTA_FILE = FASTA_FILE, SUFFIX_BAM = gsub(".bam","", SUFFIX_BAM))

saveMetricsBIN(PATH_INITIAL = PATH_INITIAL, 
               sample = SAMPLE,
               NUM_THREADS = NUM_THREADS)

# Define mapping of class to frequency and file path
CLASS_PARAMS_WGS <- list(
  Colon          = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Colorectal_Carcinoma_NCIT_C2955.tsv")),
  Lung           = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Lung_Carcinoma_NCIT_C4878.tsv")),
  Prostate       = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Prostate_Carcinoma_NCIT_C4863.tsv")),
  Urothelial     = list(freq = 22.5,file = paste0(PATH_INITIAL, "data/progenetix/Urothelial_Carcinoma_NCIT_C4030.tsv")),
  EwS            = list(freq = 10,  file = paste0(PATH_INITIAL, "data/progenetix/Ewing_Sarcoma_NCIT_C4817.tsv")),
  Mesotelioma    = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Malignant_Mesothelioma_NCIT_C4456.tsv")),
  Melanoma       = list(freq = 20,  file = paste0(PATH_INITIAL, "data/progenetix/Melanoma_NCIT_C3224.tsv")),
  Breast         = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Breast_Carcinoma_NCIT_C4872.tsv")),
  Gliobastoma    = list(freq = 35,  file = paste0(PATH_INITIAL, "data/progenetix/Gliobastoma_NCIT_C3058.tsv")),
  Astrocytoma    = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Astrocytoma_NCIT_C60781.tsv")),
  Oligodendroglioma = list(freq = 25, file = paste0(PATH_INITIAL, "data/progenetix/Oligodendroglioma_NCIT_C3288.tsv")),
  MPNST          = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/MPNST_NCIT_C3798.tsv")),
  MM             = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Multiple_Myeloma_NCIT_C3242.tsv")),
  Pancreatic     = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Pancreatic_Adenocarcinoma_NCIT_C8294.tsv"))
)

feat_WGS <- getFeatureBasedOnCNV(c(AllSample, AllSample), PATH_INITIAL = PATH_INITIAL, 
                                 CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                 NUM_THREADS = 30)


#### Feature WGS and Medip for each cancer type and healthy samples (respect to specific region) (from /home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/1_getFeat.R) FOR FIGURE 4 ####

source("lpWGS_features.R")
source("cfMeDIP_features.R")

setwd("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE")

ADAPTIVE = TRUE
WITH_MEDIP = TRUE

for(CLASS in ALL_CLASS){
  
  CLASS_CNV <- CLASS
  AllSample_OK_sub <- AllSample_OK[AllSample_OK$Class %in% c("Healthy", CLASS),]
  
  file_output <- paste0(CLASS,  "_new_feat_WGS_adaptive.RData")
  
  if(!file.exists(file_output)){
    
    wgs_mtx <- getFeatureBasedOnCNV(AllSample_OK_sub, CLASS, MIN_MAX_NORM = FALSE, PAR_CORES = 50)
    
    colnames(wgs_mtx) <- paste0( CLASS, colnames(wgs_mtx) )
    save(wgs_mtx,AllSample_OK_sub, file = file_output)
  }else{
    load(file_output)
  }
  
  
  if(WITH_MEDIP){
    
    file_output <- paste0(CLASS,  "_new_feat_WGS_adaptive_with_medip.RData")
    
    if(!file.exists(file_output)){
      medip_mtx <- getFeatMEDIP_last_COUNT(AllSample_OK_sub, MIN_MAX_NORM = FALSE)
      
      features_mtx <- cbind(wgs_mtx, medip_mtx)
      save(wgs_mtx, AllSample_OK_sub, file = file_output)
    }else{
      load(file_output)
    }
  }
  
  features_mtx
  
} 




##### CLASSIFY #####

setwd("/home3/adefalco/Fate-AI/Class_New_Fig4/")
source("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/function_adaptive.R")

#getLabels_NEW
meta <- getUpdate_metadata()

ONLY_EXTERNAL <- FALSE
ONLY_INTERNAL <- FALSE

ADAPTIVE <- TRUE

AGGREGATE_EXTERNAL_DATASET <- FALSE
#"Lung","Melanoma", "Mesotelioma","Urothelial","Colon", "Breast"
for(CLASS in c("Prostate")){
  
  for(FILTER_BIOGEM_TF in c(FALSE, TRUE)){
    
    if(CLASS %in% c("Melanoma","Mesotelioma", "Urothelial")){
      DATASETs <- c("Biogem")
      ONLY_INTERNAL <- TRUE
    }else{
      DATASETs <- c("Biogem", "GSE243474")
      ONLY_INTERNAL <- FALSE
    }
    
    if(!FILTER_BIOGEM_TF){
      ONLY_INTERNAL <- TRUE
      DATASETs <- c("Biogem")
    }
    
    if(CLASS %in% c("Prostate")){
      ONLY_INTERNAL <- TRUE
      DATASETs <- c("GSE243474")
    }
    
    AllSample_CLASSIFIER <- getSamplesWithMetaNew(CLASS, DATASETs = DATASETs, MEDIP = TRUE, FILTER_BIOGEM_TF = FILTER_BIOGEM_TF)
    
    #AllSample_CLASSIFIER[,!colSums(is.na(AllSample_CLASSIFIER))==nrow(AllSample_CLASSIFIER)]
    
    if(CLASS == "Urothelial"){
      metaURO <- as.data.frame(readxl::read_xlsx("/home3/adefalco/Fate-AI/Urothelial_Lung_filtered_samples/UrothelialClean.xlsx"))
      rownames(metaURO) <- metaURO$`ID- IC`
      
      AllSample_CLASSIFIER <- AllSample_CLASSIFIER[!rownames(AllSample_CLASSIFIER) %in% metaURO$`ID- IC`[grepl("Indenne|Assenza|Flogosi", metaURO$Diagnosi)],]
    }else if(CLASS == "Lung"){
      
      AllSample_CLASSIFIER <- AllSample_CLASSIFIER[!AllSample_CLASSIFIER$Diagnosi %in% c("TIMOMA",
                                                                                         "No Neoplasia",
                                                                                         "Non diagnostico",
                                                                                         "TIMOMA",
                                                                                         "sarcoma mixoide"
                                                                                         #,"No Procedura"
      ),]
      
    }
    
    print("TABLE CLASS DATASET")
    print(table(AllSample_CLASSIFIER$Class, AllSample_CLASSIFIER$Dataset))
    
    if(ONLY_INTERNAL & FILTER_BIOGEM_TF) next
    
    #AllSample_CLASSIFIER <- AllSample_CLASSIFIER[!AllSample_CLASSIFIER$Diagnosi %in% c("TIMOMA","sarcoma mixoide", "Non diagnostico"),]
    
    for(METODO in c("delfiMEDIP","FateAI", "FateAI_medip", "GriffinMEDIP")){   #"FateAI_medip", "FateAI", "delfiMEDIP", "onlyMEDIP_NEW", GriffinMEDIP")){       
      
      print(METODO)
      
      if(METODO=="GriffinMEDIP" & FILTER_BIOGEM_TF) next
      
      if(METODO=="GriffinMEDIP" & Class=="Prostate") next
      
      AllSample <- AllSample_CLASSIFIER
      
      if(ADAPTIVE){
        load(paste0("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE/",CLASS,"_new_feat_WGS_adaptive_with_medip.RData"))
      }else{
        load(paste0("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE/",CLASS,"_new_feat_WGS_with_medip.RData"))
      }
      
      #load("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE/Lung_new_feat_WGS_adaptive.RData")
      AllSample_OK <- AllSample_OK_sub
      rownames(AllSample_OK) <- AllSample_OK$IC_code
      rownames(delfi_mtx) <- rownames(AllSample_OK)
      delfi_mtx <- delfi_mtx[rownames(AllSample),]
      medip_mtx <- delfi_mtx[, grepl("Ratio|Score", colnames(delfi_mtx))]
      delfi_mtx <- delfi_mtx[, !grepl("Ratio|Score", colnames(delfi_mtx))]
      
      delfix_mtx <- delfi_mtx[, grepl(CLASS, colnames(delfi_mtx))]
      
      
      if(METODO %in% c("delfi", "delfiMEDIP", "delfiMEDIP_cor")){
        sampNAME <- AllSample$FASTQ_Name
        sampNAME[grepl("ICA", sampNAME)] <- gsub("_S[0-9][0-9]","", sampNAME[grepl("ICA", sampNAME)] )
        sampNAME[grepl("ICA", sampNAME)] <- gsub("_S[0-9]","", sampNAME[grepl("ICA", sampNAME)] )
        
        sampNAME[grepl("Sample", sampNAME)] <- AllSample[grepl("Sample", AllSample$FASTQ_Name),]$IC_code 
        
        sampNAME[grepl("Grisolia-18731", sampNAME)] <- AllSample[grepl("Grisolia-18731", AllSample$FASTQ_Name),]$IC_code 
        
        sampNAME[AllSample$Dataset=="GSE243474"] <- AllSample[AllSample$Dataset=="GSE243474",]$IC_code
        
        delfi_mtx <- getDelfiFeatures(sampNAME)
        
        if(length(setdiff(sampNAME, rownames(delfi_mtx)))>0) stop("error message")
        
        delfi_mtx <- delfi_mtx[sampNAME,]
        rownames(delfi_mtx) <- AllSample$FASTQ_Name
        
        METHOD_CLASSIFIER <- "gbm"
      }else if(METODO == "FateAI"){
        
        
        delfi_mtx <- delfi_mtx
        
        #delfi_mtx <- getFeatureBasedOnCN(AllSample, CLASS, MIN_MAX_NORM = TRUE)
        METHOD_CLASSIFIER <- "glmnet"
        
        delfi_mtx <- normalizeMatrixDataset(delfi_mtx, AllSample)
        
      }else if(METODO == "FateAI_medip"){
        
        #delfi_mtx <- getFeatureBasedOnCN(AllSample, CLASS, MIN_MAX_NORM = TRUE)
        #medip_mtx <- getFeatMEDIP_last_COUNT(AllSample, MIN_MAX_NORM = TRUE)
        delfi_mtx <- cbind(delfi_mtx, medip_mtx)
        
        METHOD_CLASSIFIER <- "glmnet"
        
        delfi_mtx <- normalizeMatrixDataset(delfi_mtx, AllSample)
        
      }else if (METODO == "FATE_MEDIP_TEST"){
        
        delfi_mtx <- cbind(delfi_mtx, medip_mtx)
        
        delfi_mtx <- normalizeMatrixDataset(delfi_mtx, AllSample)
        
        METHOD_CLASSIFIER <- "svmLinear"
        #svmLinear svmRadial svmPoly
      }else if(METODO == "onlyMEDIP_NEW"){
        
        #medip_mtx <- getFeatMEDIP_last_COUNT(AllSample, MIN_MAX_NORM = TRUE)
        delfi_mtx <- medip_mtx
        
        METHOD_CLASSIFIER <- "glmnet"  
        
        delfi_mtx <- normalizeMatrixDataset(delfi_mtx, AllSample)
        
      }else if (METODO %in% c("Griffin", "GriffinMEDIP")){
        
        sampNAME <- AllSample$FASTQ_Name
        sampNAME[grepl("ICA", sampNAME)] <- gsub("_S[0-9][0-9]","", sampNAME[grepl("ICA", sampNAME)] )
        sampNAME[grepl("ICA", sampNAME)] <- gsub("_S[0-9]","", sampNAME[grepl("ICA", sampNAME)] )
        
        sampNAME[grepl("Sample", sampNAME)] <- AllSample[grepl("Sample", AllSample$FASTQ_Name),]$IC_code 
        
        sampNAME[grepl("Grisolia-18731", sampNAME)] <- AllSample[grepl("Grisolia-18731", AllSample$FASTQ_Name),]$IC_code 
        
        delfi_mtx <- getGriffinFeatures(sampNAME)
        
        if(length(setdiff(sampNAME, rownames(delfi_mtx)))>0) stop("error message")
        
        delfi_mtx <- delfi_mtx[sampNAME,]
        rownames(delfi_mtx) <- AllSample$FASTQ_Name
        
        #delfi_mtx <- getGriffinFeatures(AllSample$FASTQ_Name)
        delfi_mtx$type <- NULL
        VAR_LIMIT = 0.80
        resProf_train_pca <- prcomp(delfi_mtx, center = TRUE, scale. = TRUE)
        res_pca <- summary(resProf_train_pca)
        res_pca_var <- cumsum(res_pca$importance[2,])
        numPC <- which.min(res_pca_var <= VAR_LIMIT)
        #print(res_pca_var[numPC])
        delfi_mtx <- resProf_train_pca$x[,1:numPC]
        
        METHOD_CLASSIFIER <- "glmnet"
      }
      
      file_name = paste(METODO, CLASS, METHOD_CLASSIFIER, sep = "_")
      if(FILTER_BIOGEM_TF) file_name <- paste0(file_name, "_tf_filtered")
      file_name <- paste0(file_name, ".RData")
      
      setwd("/home3/adefalco/Fate-AI/Class_New_Fig4/SAVE_FILE_CLASSIFIER/")
      
      
      if(!file.exists(file_name)){
        
        MTX <- as.data.frame(delfi_mtx)
        
        pheatmap::pheatmap(MTX)
        
        combinations <- expand.grid(DATASETs,DATASETs)
        colnames(combinations) <- c("TRAIN", "TEST")
        
        if(ONLY_EXTERNAL) {
          combinations <- combinations[!apply(combinations, 1, function(x) x[1]==x[2]),]
        }
        
        
        
        if(ONLY_INTERNAL) {
          combinations <- combinations[apply(combinations, 1, function(x) x[1]==x[2]),]
        }
        
        
        
        if(FILTER_BIOGEM_TF){
          combinations <- combinations[!combinations$TEST=="Biogem",]
        }
        
        
        
        if(!FILTER_BIOGEM_TF & ONLY_EXTERNAL){
          combinations <- combinations[!combinations$TRAIN=="Biogem",]
          combinations <- combinations[combinations$TEST=="Biogem",]
        }
        
        #MTX$type <- factor(MTX$type, levels = c(unique(MTX$type)[unique(MTX$type)!="Healthy"],unique(MTX$type)[unique(MTX$type)=="Healthy"]))
        
        getAUC <- function(resPred, label0){
          library(ROCR)
          labels <- ifelse(resPred[,ncol(resPred)] == label0, 1, 0)
          scores <- resPred[,1]
          pred <- prediction(scores, labels)
          perf <- performance(pred, "tpr", "fpr")
          
          AUC <- performance(pred, "auc")
          
          AUC <- AUC@y.values[[1]]
          AUC
        }
        
        
        if(class(combinations)!="data.frame"){
          num_elem <- length(combinations)
        }else{
          num_elem <- nrow(combinations)
        }
        
        print(paste0("combinations"))
        print(combinations)
        
        resAUCs <- parallel::mclapply(1:num_elem, function(i){
          
          if(class(combinations)!="data.frame"){
            combinations_sub <- combinations[[i]]
          }else{
            combinations_sub <- combinations[i,]
          }
          
          
          print(combinations_sub)
          MTX_sub <- MTX[which(AllSample$Dataset %in% c(combinations_sub$TRAIN, combinations_sub$TEST)),]
          AllSample_SUB <- AllSample[which(AllSample$Dataset %in% c(combinations_sub$TRAIN, combinations_sub$TEST)),] 
          
          if(any(combinations_sub$TRAIN %in% combinations_sub$TEST)){
            
            res_delfi <- classifyMATRIX(MTX_sub, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, method = METHOD_CLASSIFIER)
            TYPE = "100x10fold-CV"
            
            TRAIN_SAMP <- paste(rownames(AllSample_SUB), collapse = ";")
            TEST_SAMP <- paste(rownames(AllSample_SUB), collapse = ";")
            NUM_CLASS_TRAIN <- table(AllSample_SUB$Class)
            NUM_CLASS_TEST <- table(AllSample_SUB$Class)
            
            NUM_CLASS_TRAIN <- paste(paste(names(NUM_CLASS_TRAIN), NUM_CLASS_TRAIN, sep = ":"), collapse = ";")
            NUM_CLASS_TEST <- paste(paste(names(NUM_CLASS_TEST), NUM_CLASS_TEST, sep = ":"), collapse = ";")
            
          }else{
            
            TEST_INDEX <- which(AllSample_SUB$Dataset %in% combinations_sub$TEST)
            res_delfi <- classifyMATRIX(MTX_sub, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, testIND = TEST_INDEX, method = METHOD_CLASSIFIER)
            TYPE = "EXTERNAL_VALIDATION"
            
            TRAIN_SAMP <- paste(rownames(AllSample_SUB[setdiff(1:nrow(AllSample_SUB),TEST_INDEX),]), collapse = ";")
            TEST_SAMP <- paste(rownames(AllSample_SUB[TEST_INDEX,]), collapse = ";")
            NUM_CLASS_TRAIN <- table(AllSample_SUB[setdiff(1:nrow(AllSample_SUB),TEST_INDEX),]$Class)
            NUM_CLASS_TEST <- table(AllSample_SUB[TEST_INDEX,]$Class)
            
            NUM_CLASS_TRAIN <- paste(paste(names(NUM_CLASS_TRAIN), NUM_CLASS_TRAIN, sep = ":"), collapse = ";")
            NUM_CLASS_TEST <- paste(paste(names(NUM_CLASS_TEST), NUM_CLASS_TEST, sep = ":"), collapse = ";")
            
          }
          
          AUC_gbm <- getAUC(res_delfi, label0 = CLASS)
          
          if(length(combinations_sub$TEST)>1){
            nameTEST <- paste(combinations_sub$TEST, sep = "_")
          }else{
            nameTEST <- combinations_sub$TEST
          }
          
          if(length(combinations_sub$TRAIN)>1){
            nameTRAIN <- paste(combinations_sub$TRAIN, sep = "_")
          }else{
            nameTRAIN <- combinations_sub$TRAIN
          }
          
          resAUCs <- data.frame(TRAIN = nameTRAIN, TEST = nameTEST, gbm = AUC_gbm, TYPE = TYPE, METHOD_CLASSIFIER = METHOD_CLASSIFIER, METODO = METODO)
          
          resAUCs <- cbind(resAUCs, TRAIN_SAMP, TEST_SAMP, NUM_CLASS_TRAIN, NUM_CLASS_TEST, FILTER_BIOGEM_TF)
          
          #resAUCs
          list(res_delfi, resAUCs)
        }, mc.cores = 1)
        
        
        #if(AGGREGATE_EXTERNAL_DATASET) CLASS <- paste0(CLASS, "_AGGREGATE_DATASETS")
        

        save(resAUCs,MTX, AllSample,combinations, file = file_name)
        
        
      }
      
    }
    
  }
}

