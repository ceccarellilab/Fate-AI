
# Analysis of DNA methylation data
# Analysis of TCGA 450K methylation data was done with TCGABioloinks (v.2.12.6)62. 
# Data for normal colon and colorectal adenocarcinoma were downloaded using the function GDCquery with project=c(‘TCGA-COAD’), 
# data.category=‘DNA Methylation’, legacy=FALSE, platform=c(‘Illumina Human Methylation 450’) and 
# sample.type=c(‘Primary solid Tumor’,’Solid Tissue Normal’). Differentially methylated probes between 
# normal and colorectal adenocarcinoma samples were computed using TCGAanalyze_DMR with a P value cutoff of 10−5 and a 
# mean difference in β value cutoff of 0.25 to determine significance. Overlaps between DNA methylation probes and our 
# peak set were identified with the GenomicRanges function FindOverlaps in R.

#### DMRs analysis from TCGA ####

saveDMRs_fromTCGA <- function(PATH_INITIAL = "./", CancerTypes = c("TCGA-BRCA", "TCGA-LUAD","TCGA-LUSC","TCGA-PRAD", "TCGA-BLCA", "TCGA-SKCM", "TCGA-MESO"), NUM_THREADS = 1, TCGA_DIR = "data/tcga/"){

  library(parallel)

  parallel::mclapply(CancerTypes, function(TCGA_CLASS){
  
    WORK_DIR <- paste0(PATH_INITIAL, TCGA_DIR, TCGA_CLASS)
    
    METH_PATH <- paste0(WORK_DIR, "/", TCGA_CLASS, "_METH.RData")
    DMR_PATH <- paste0(PATH_INITIAL, TCGA_DIR, "/", TCGA_CLASS, "_METH_data_DMR.RData")
    
    if (!dir.exists(WORK_DIR)) {
      dir.create(WORK_DIR)
    }
    setwd(WORK_DIR)
    
    library(TCGAbiolinks)
    
    if(!file.exists(METH_PATH) & !file.exists(DMR_PATH)){
      
      message("Download Meth Data")
      
      query <- GDCquery(
        project = TCGA_CLASS,
        data.category = "DNA Methylation",
        platform = "Illumina Human Methylation 450",
        sample.type = c("Primary Tumor","Solid Tissue Normal"),
        data.type = "Methylation Beta Value"
      )
      
      table(query$results[[1]]$sample_type)
      
      GDCdownload(query, method = "api", files.per.chunk = 40)
      data <- GDCprepare(query)
      
      save(data, file = METH_PATH)
    }
    
    if(!file.exists(DMR_PATH)){
      
      message("DMR Analysis")
      
      load(METH_PATH)
    
      library(SummarizedExperiment)
      
      data.met <- data
      data.met <- subset(data.met,subset = (rowSums(is.na(assay(data.met))) == 0))
      
      data_DMR <- TCGAanalyze_DMC(
        data = data.met, 
        groupCol = "sample_type",
        group1 = "Primary Tumor",
        group2 = "Solid Tissue Normal",
        p.cut = 10^-5,
        diffmean.cut = 0.25,
        legend = "State",
        plot.filename = paste0(TCGA_CLASS,"_volcano.png")
      )
      
      save(data_DMR, file = DMR_PATH)
    
    }
    
    }, mc.cores = NUM_THREADS)
  
}

######## Select top hypermethylated from TCGA and MethAtlas ########

saveBED_TopDMRs <- function(PATH_INITIAL = "./", ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma"), N_TOP = 3000, sizeBIN = 300, TCGA_DIR = "data/tcga/", BED_DIR = "data/bed_dmrs/"){

  #Forward = Plus = (+) = Positive = Watson
  #Reverse = Minus = (-) = Negative = Crick
  
  library(deconvR) 
  
  data("HumanCellTypeMethAtlas")
  
  data("IlluminaMethEpicB5ProbeIDs")
  
  ALTs <- paste0("Hyper_", ClassTypes)
  
  data_DMR <- lapply(ALTs, function(CLASS){
    
    if(CLASS == "Hyper_Plasma"){
      Healthy_Types <- c("Vascular_endothelial_cells","Hepatocytes","Erythrocyte_progenitors","Monocytes_EPIC","Neutrophils_EPIC","B.cells_EPIC","CD4T.cells_EPIC","CD8T.cells_EPIC","NK.cells_EPIC")
      rownames(HumanCellTypeMethAtlas) <- HumanCellTypeMethAtlas$IDs
      HumanCellTypeMethAtlas <- HumanCellTypeMethAtlas[,-1]
      
      Score <- apply(HumanCellTypeMethAtlas[,colnames(HumanCellTypeMethAtlas) %in% Healthy_Types], 1, mean) - apply(HumanCellTypeMethAtlas[,!colnames(HumanCellTypeMethAtlas) %in% Healthy_Types], 1, mean)
      
      Score <- Score[order(Score, decreasing = TRUE)]
      Score <- Score[Score>0]
      top_Plasma <- names(Score)[1:min(length(Score), N_TOP)]
      top_Plasma
    }else{
      
      if (!CLASS %in% names(CLASS_TO_TCGA)) {
        stop(paste("Unknown CLASS:", CLASS))
      }
      
      tcga_code <- CLASS_TO_TCGA[[CLASS]]
      file_path <- paste0(PATH_INITIAL, TCGA_DIR, tcga_code, "_METH_data_DMR.RData")
      
      # Caricamento file
      load(file_path)
      
      nameHYPER <- "Hypermethylated in Primary Tumor"
      data_DMR <- data_DMR[data_DMR$status==nameHYPER,]
      data_DMR$FC <- log2(data_DMR$mean.Primary.Tumor/data_DMR$mean.Solid.Tissue.Normal)
      data_DMR <- data_DMR[data_DMR$FC>1.5 & data_DMR$p.value.adj.Primary.Tumor.Solid.Tissue.Normal<0.01,]
      data_DMR$SCORE <- -log10(data_DMR$p.value.adj.Primary.Tumor.Solid.Tissue.Normal)*data_DMR$FC
      data_DMR <- data_DMR[order(data_DMR$SCORE, decreasing = TRUE),]
      
      data_DMR <- data_DMR[1:min(nrow(data_DMR), N_TOP),]
      top_Tumor <- rownames(data_DMR)
      top_Tumor
    }
    
  })
  
  names(data_DMR) <-   ALTs 
  
  Hyper_Lung_SHARED <- intersect(data_DMR[["Hyper_Lung_LUAD"]], data_DMR[["Hyper_Lung_LUSC"]])
  
  tabCOUNT <- table(unlist(data_DMR[-which(names(data_DMR)=="Hyper_Plasma")]))
  Hyper_TUMOR_SHARED <- names(tabCOUNT[tabCOUNT>2])
  
  Hyper_Lung_SHARED <- Hyper_Lung_SHARED[!Hyper_Lung_SHARED %in% Hyper_TUMOR_SHARED]
  
  data_DMR <- lapply(1:length(data_DMR), function(n) setdiff(data_DMR[[n]], unlist(data_DMR[-n])))    
  names(data_DMR) <-   ALTs 
  
  lapply(data_DMR, function(x) length(x))
  
  data_DMR$Hyper_Lung_SHARED <- Hyper_Lung_SHARED
  data_DMR$Hyper_TUMOR_SHARED <- Hyper_TUMOR_SHARED
  ALTs <- names(data_DMR)
  
  annotation <- read.table(paste0(PATH_INITIAL,"data/HM450.hg38.manifest.tsv.gz"), header = TRUE, sep = "\t")
  
  resFEATUREs_alt <- lapply(ALTs, function(ALT){
    
    PATH_WATSON <- paste0(PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Watson.bed")
    PATH_CRICK <- paste0(PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Crick.bed")
    
    if(!file.exists(PATH_WATSON) | !file.exists(PATH_CRICK)){
      df_BED <- annotation[annotation$probeID %in% (data_DMR[[ALT]]),1:4]
      centerBIN <- df_BED$CpG_end-1
      df_BED$CpG_beg <- centerBIN - sizeBIN/2
      df_BED$CpG_end <- centerBIN + sizeBIN/2
      ALT <- gsub(" in Primary Tumor", "",ALT )
      
      write.table(df_BED[df_BED$probe_strand=="+",1:3], file = PATH_WATSON, quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(df_BED[df_BED$probe_strand=="-",1:3], file = PATH_CRICK, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  
  })

}


######### Extract coverage in DMRs from cfMedip bam #########

saveCoverageDMRs_fromBam <- function(PATH_INITIAL = "./", ALL_BAM_MEDIP_PATH = "/home3/adefalco/Fate-AI/ScriptMedip/BAM_MEDIP/", 
                                     FASTA_FILE = "/storage/qnap_vol1/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa",
                                     PATH_SAMTOOLS = "/home/adefalco/singleCell/cellRank/samtools-1.11/samtools",
                                     ClassTypes = c("Colon", "Lung_LUAD", "Lung_LUSC","Breast", "Prostate","Urothelial","Melanoma", "Mesotelioma", "Plasma", "Lung_SHARED"),
                                     SUFFIX_BAM = ".sorted.bam", 
                                     N_TOP = 3000, 
                                     NUM_THREADS = 10,   
                                     BIN_SIZE = 300,
                                     MAPQ = 30,
                                     MAX_FRAG_LENGHT = 550,
                                     MIN_FRAG = 20,
                                     SUFFIX_SAVE_FILE = "_res_frag_motif.RData",
                                     tmp_dir = "/temp",
                                     awk_file_filter = "scripts/filter_cfmedip.awk",
                                     BED_DIR = "data/bed_dmrs/",
                                     DMR_COUNT_DIR = "output/cfMeDIP/"){

  library(GenomicRanges)

  dirSave = paste0(PATH_INITIAL, DMR_COUNT_DIR, "FRAGM_EXTR_",N_TOP, "_COUNTS/")

  ALL_BAM_MEDIP <- getPathBam(MEDIP = T)
  AllSample <- getSamples(MEDIP = T)
  
  ALTs <- paste0("Hyper_", ClassTypes)
  
  COUNTS_samples <- parallel::mclapply(1:length(AllSample), function(i){
    
    sample <- AllSample[i]
    bam <- ALL_BAM_MEDIP[i]
      
      path_output <- paste0(dirSave,sample,"_",as.integer(BIN_SIZE),SUFFIX_SAVE_FILE)
      
      if(!file.exists(path_output)){
        
        dir.create(dirSave)
        
        library(parallel)
        library(dplyr)
        
        #setwd("../")
        
        #setwd("../")
        #print(getwd())
        
        resFEATUREs_alt <- lapply(ALTs, function(ALT){
          
          watson <- paste0(PATH_SAMTOOLS, " view ", bam, " -f 99 ", "-M -L ", PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Watson.bed", " | awk -v MIN_MAPQ=", MAPQ, " -v MAX_FRAGMENT_LEN=", MAX_FRAG_LENGHT, " -f ", paste0(PATH_INITIAL,awk_file_filter))
          
          crick <- paste0(PATH_SAMTOOLS, " view ", bam, " -f 163 ", "-M -L ", PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Crick.bed", " | awk -v MIN_MAPQ=", MAPQ, " -v MAX_FRAGMENT_LEN=", MAX_FRAG_LENGHT, " -f ", paste0(PATH_INITIAL,awk_file_filter))
          
          #### WATSON
          
          watsonFRAG <- tryCatch( 
            {
              read.csv(text = system(watson, intern = TRUE), 
                       header = FALSE, sep = "\t")
            },
            error = function(e) {
              NULL
            }
          )
          
          crickFRAG <- tryCatch( 
            {
              read.csv(text = system(crick, intern = TRUE), 
                       header = FALSE, sep = "\t")
            },
            error = function(e) {
              NULL
            }
          )
          
          
          if (!is.null(watsonFRAG)){
            watsonFRAG = watsonFRAG[watsonFRAG$V9>=MIN_FRAG,]
            watsonFRAG$V1 <- "watsonFRAG"
          }
          if (!is.null(crickFRAG)){
            crickFRAG = crickFRAG[crickFRAG$V9>=MIN_FRAG,]
            crickFRAG$V1 <- "crickFRAG"
          }
          
          watsonFRAG_GR <- GenomicRanges::makeGRangesFromDataFrame(watsonFRAG, seqnames.field = "V2", start.field = "V3", end.field = "V4")
          
          bed_GR <- GenomicRanges::makeGRangesFromDataFrame(read.table(paste0(PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Watson.bed")), seqnames.field = "V1", start.field = "V2", end.field = "V3")
          bed_GR <- reduce(resize(bed_GR, width(bed_GR), "start"))
          
          counts_watson <- table(findOverlaps(watsonFRAG_GR, bed_GR)@to)
          names(counts_watson) <- paste0(names(counts_watson),"_", ALT, "_TOP", N_TOP, "_Watson.bed")
          
          crickFRAG_GR <- GenomicRanges::makeGRangesFromDataFrame(crickFRAG, seqnames.field = "V2", start.field = "V3", end.field = "V4")
          
          bed_GR <- GenomicRanges::makeGRangesFromDataFrame(read.table(paste0(PATH_INITIAL, BED_DIR,ALT, "_TOP", N_TOP, "_Crick.bed")), seqnames.field = "V1", start.field = "V2", end.field = "V3")
          bed_GR <- reduce(resize(bed_GR, width(bed_GR), "start"))
          
          counts_crick <- table(findOverlaps(crickFRAG_GR, bed_GR)@to)
          names(counts_crick) <- paste0(names(counts_crick),"_", ALT, "_TOP", N_TOP, "_Crick.bed")
          
          counts <- c(counts_watson, counts_crick)
          counts <- as.data.frame(counts)
          colnames(counts) <- sample
          counts
        })
        
        resFEATUREs_alt <- Reduce(rbind, resFEATUREs_alt)
        
        save(resFEATUREs_alt, file = path_output)
      }else{
        load(path_output)
      }
    
    resFEATUREs_alt
    
  }, mc.cores = NUM_THREADS) 

}

##### Aggregate coverage of all samples (save) ####

merge_Samples_cfMEDIP_Counts <- function(PATH_INITIAL = "./" , 
                                         N_TOP = 3000,
                                         DMR_COUNT_DIR = "output/cfMeDIP/"){
  
  dirLoad = paste0(PATH_INITIAL, DMR_COUNT_DIR, "FRAGM_EXTR_",N_TOP, "_COUNTS/")
  DMR_COUNT_MERGE_PATH <-  paste0(PATH_INITIAL, DMR_COUNT_DIR, "COUNTS_samples_merge_MEDIP_", N_TOP, ".RData")
  
  Res_medip <- lapply(list.files(dirLoad), function(x) get(load(paste0(dirLoad,x))))
  
  COUNTS_samples <- lapply(Res_medip, function(COUNT){
    COUNT$ID <- rownames(COUNT)
    COUNT
  })
  
  COUNTS_samples_merge <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all = TRUE),
                                 COUNTS_samples)
  COUNTS_samples_merge[is.na(COUNTS_samples_merge)] <- 0
  
  rownames(COUNTS_samples_merge) <- COUNTS_samples_merge$ID
  COUNTS_samples_merge$ID <- NULL
  
  save(COUNTS_samples_merge, file = DMR_COUNT_MERGE_PATH)
}

##### Get features cfMedip  ####

getFeature_cfMeDIP <- function(AllSample, PATH_INITIAL = "./", ALL_BAM_MEDIP_PATH = "/home3/adefalco/Fate-AI/ScriptMedip/BAM_MEDIP/", SUFFIX_BAM = ".sorted.bam",  ClassTypes = c("Plasma", "Colon", "Prostate", "Breast", "Lung", "Mesotelioma", "Melanoma", "Urothelial"), 
                               N_TOP = 3000,
                               DMR_COUNT_DIR = "output/cfMeDIP/"){
  
  DMR_COUNT_MERGE_PATH <-  paste0(PATH_INITIAL, DMR_COUNT_DIR, "COUNTS_samples_merge_MEDIP_", N_TOP, ".RData")
  
  load(DMR_COUNT_MERGE_PATH)
  
  COUNTS_samples_merge[COUNTS_samples_merge<3] <- 0
  
  score_classes <- lapply(ClassTypes, function(CLASS){
    SUBSET_MTX <- COUNTS_samples_merge[grepl(CLASS, rownames(COUNTS_samples_merge)),]
    score <- apply(SUBSET_MTX, 2, function(x) sum(x>0))/nrow(SUBSET_MTX)
  })
  
  score_classes <- as.data.frame(Reduce(cbind, score_classes))
  colnames(score_classes) <- paste0("Score", ClassTypes)
  
  coverage_classes <- lapply(ClassTypes, function(CLASS){
    SUBSET_MTX <- COUNTS_samples_merge[grepl(CLASS, rownames(COUNTS_samples_merge)),]
    score <- apply(SUBSET_MTX, 2, function(x) sum(x))/nrow(SUBSET_MTX)
  })
  
  names(coverage_classes) <- ClassTypes
  
  ratioCoverageList <- lapply(coverage_classes[-1], function(score){
    score/coverage_classes[[1]]
  })
  ratioCoverage <- as.data.frame(Reduce(cbind, ratioCoverageList))
  colnames(ratioCoverage) <- paste0("Ratio", names(ratioCoverageList))
  
  feat_medip <- cbind(score_classes, ratioCoverage)

  feat_medip[AllSample,]
  
}
