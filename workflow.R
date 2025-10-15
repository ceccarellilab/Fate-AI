##### WORKFLOW cfMEDIP#####

FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools"
NUM_THREADS <- 10
PATH_INITIAL <- "/home2/adefalco/Fate-AI"
ALL_BAM_MEDIP_PATH = "/home3/adefalco/Fate-AI/FIGURE/FIGURES_paper/test_script/TEST_BAM_MEDIP"
SUFFIX_BAM = ".sorted.bam"

#TCGA_DIR <- "data/tcga/"
#BED_DIR <- "data/bed_dmrs"
#N_TOP <- 3000

source(paste0(PATH_INITIAL, "scripts/cfMeDIP_features.R"))

CLASS_TO_TCGA <- list(
  "Hyper_Colon"       = "TCGA_COAD",
  "Hyper_Lung_LUAD"   = "TCGA-LUAD",
  "Hyper_Lung_LUSC"   = "TCGA-LUSC",
  "Hyper_Breast"      = "TCGA-BRCA",
  "Hyper_Prostate"    = "TCGA-PRAD",
  "Hyper_Melanoma"    = "TCGA-SKCM",
  "Hyper_Mesotelioma" = "TCGA-MESO",
  "Hyper_Urothelial"  = "TCGA-BLCA"
)


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
  MM             = list(freq = 25,  file = paste0(PATH_INITIAL, "data/progenetix/Multiple_Myeloma_NCIT_C3242.tsv")),
  Pancreatic     = list(freq = 15,  file = paste0(PATH_INITIAL, "data/progenetix/Pancreatic_Adenocarcinoma_NCIT_C8294.tsv"))
)

feat_WGS <- getFeatureBasedOnCNV(AllSample, PATH_INITIAL = PATH_INITIAL, 
                                 CLASS_CNV = names(CLASS_PARAMS_WGS)[1], 
                                 NUM_THREADS = 30)


#### Feature WGS and Medip ####

source(paste0(PATH_INITIAL, "scripts/cfMeDIP_features.R"))
source(paste0(PATH_INITIAL, "scripts/lpWGS_features.R"))

CLASS <- "Colon"
METHOD_CLASSIFIER <- "glmnet"
MODEL <- "Fate-AI(+Meth)"

      load(paste0("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE/",CLASS,"_new_feat_WGS_adaptive_with_medip.RData"))
 
      #load("/home3/adefalco/Fate-AI/ADAPTIVE_ANALYSIS/SAVE_FILE/Lung_new_feat_WGS_adaptive.RData")
      AllSample_OK <- AllSample_OK_sub
      rownames(AllSample_OK) <- AllSample_OK$IC_code
      rownames(feat_mtx) <- rownames(AllSample_OK)
      feat_mtx <- feat_mtx[rownames(AllSample),]
      medip_mtx <- feat_mtx[, grepl("Ratio|Score", colnames(feat_mtx))]
      feat_mtx <- feat_mtx[, !grepl("Ratio|Score", colnames(feat_mtx))]
      
      delfix_mtx <- feat_mtx[, grepl(CLASS, colnames(feat_mtx))]
      
      if(MODEL == "Fate-AI"){
        
        feat_mtx <- getFeatureBasedOnCN(AllSample, CLASS)
        METHOD_CLASSIFIER <- "glmnet"
        
        feat_mtx <- normalizeMatrixDataset(feat_mtx, AllSample)
        
      }else if(MODEL == "Fate-AI(+Meth)"){
        
        feat_mtx <- getFeatureBasedOnCN(AllSample, CLASS)
        medip_mtx <- getFeatMEDIP_last_COUNT(AllSample)
        feat_mtx <- cbind(feat_mtx, medip_mtx)
        
        METHOD_CLASSIFIER <- "glmnet"
        
        feat_mtx <- normalizeMatrixDataset(feat_mtx, AllSample)
      }
      
      file_name = paste(MODEL, CLASS, METHOD_CLASSIFIER, sep = "_")
      if(FILTER_BIOGEM_TF) file_name <- paste0(file_name, "_tf_filtered")
      file_name <- paste0(file_name, ".RData")
      
      
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
            
            prediction <- classifyMATRIX(MTX_sub, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, method = METHOD_CLASSIFIER)
            TYPE = "100x10fold-CV"
            
            TRAIN_SAMP <- paste(rownames(AllSample_SUB), collapse = ";")
            TEST_SAMP <- paste(rownames(AllSample_SUB), collapse = ";")
            NUM_CLASS_TRAIN <- table(AllSample_SUB$Class)
            NUM_CLASS_TEST <- table(AllSample_SUB$Class)
            
            NUM_CLASS_TRAIN <- paste(paste(names(NUM_CLASS_TRAIN), NUM_CLASS_TRAIN, sep = ":"), collapse = ";")
            NUM_CLASS_TEST <- paste(paste(names(NUM_CLASS_TEST), NUM_CLASS_TEST, sep = ":"), collapse = ";")
            
          }else{
            
            TEST_INDEX <- which(AllSample_SUB$Dataset %in% combinations_sub$TEST)
            prediction <- classifyMATRIX(MTX_sub, classes = AllSample_SUB$Class, class1 = "Healthy", class2 = CLASS, testIND = TEST_INDEX, method = METHOD_CLASSIFIER)
            TYPE = "EXTERNAL_VALIDATION"
            
            TRAIN_SAMP <- paste(rownames(AllSample_SUB[setdiff(1:nrow(AllSample_SUB),TEST_INDEX),]), collapse = ";")
            TEST_SAMP <- paste(rownames(AllSample_SUB[TEST_INDEX,]), collapse = ";")
            NUM_CLASS_TRAIN <- table(AllSample_SUB[setdiff(1:nrow(AllSample_SUB),TEST_INDEX),]$Class)
            NUM_CLASS_TEST <- table(AllSample_SUB[TEST_INDEX,]$Class)
            
            NUM_CLASS_TRAIN <- paste(paste(names(NUM_CLASS_TRAIN), NUM_CLASS_TRAIN, sep = ":"), collapse = ";")
            NUM_CLASS_TEST <- paste(paste(names(NUM_CLASS_TEST), NUM_CLASS_TEST, sep = ":"), collapse = ";")
            
          }
          
          AUC_gbm <- getAUC(prediction, label0 = CLASS)
          
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
          
          resAUCs <- data.frame(TRAIN = nameTRAIN, TEST = nameTEST, gbm = AUC_gbm, TYPE = TYPE, METHOD_CLASSIFIER = METHOD_CLASSIFIER, MODEL = MODEL)
          
          resAUCs <- cbind(resAUCs, TRAIN_SAMP, TEST_SAMP, NUM_CLASS_TRAIN, NUM_CLASS_TEST, FILTER_BIOGEM_TF)
          
          list(prediction, resAUCs)
        }, mc.cores = 1)
        
        
        save(resAUCs,MTX, AllSample,combinations, file = file_name)
        
        
      }
      
    }
    


