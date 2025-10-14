#Extract CNV regions from Progenetix
getCNV_Regions <- function(CLASS, FREQ_MANUAL = NULL, FREQ_MANUAL_GAIN = NULL, FREQ_MANUAL_LOSS = NULL){
  
  if (!CLASS %in% names(CLASS_PARAMS_WGS)) {
    stop("Unknown CLASS: ", CLASS)
  }
  
  params <- CLASS_PARAMS_WGS[[CLASS]]
  
  FREQ <- params$freq
  CNV_FREQ <- read.table(params$file, header = TRUE)
  
  if(!is.null(FREQ_MANUAL)){
    FREQ_GAIN <- FREQ_LOSS <- FREQ_MANUAL
  }else{
    FREQ_GAIN <- FREQ_LOSS <- FREQ
  }

  if(!is.null(FREQ_MANUAL_GAIN)) FREQ_GAIN <- FREQ_MANUAL_GAIN
  if(!is.null(FREQ_MANUAL_LOSS)) FREQ_LOSS <- FREQ_MANUAL_LOSS
  
  CNV_FREQ_ALL <- CNV_FREQ
  
  GAIN <- CNV_FREQ[CNV_FREQ$gain_frequency>=FREQ_GAIN,]
  GAIN$ALT <- "GAIN"
  
  LOSS <- CNV_FREQ[CNV_FREQ$loss_frequency>=FREQ_LOSS,]
  LOSS$ALT <- "LOSS"
  
  CNV_FREQ <- rbind(GAIN, LOSS)
  
  CNV_GR <- GenomicRanges::makeGRangesFromDataFrame(CNV_FREQ, seqnames.field = "reference_name", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

  CNV_GR <- unlist(GenomicRanges::reduce(split(CNV_GR, ~ALT)))
  CNV_GR$ALT <- names(CNV_GR)
  CNV_GR <- CNV_GR[CNV_GR@seqnames %in% 1:22,]
  
  CNV_GR
}

#Extract WGS features based on CNV regions from Progenetix
getFeatureBasedOnCNV <- function(AllSample, 
          CLASS_CNV, 
          NUM_THREADS = 30,
          MIN_MAX_NORM = FALSE, 
          FREQ = NULL, 
          MIN_SIZE_ALT = 1500000, 
          features_sel = c("mean", "ratio_NucCor_Nuc" , "ratio_NucCorChrom_Nuc", "coverage","coverageNucCore", "coverageChrom","coverageNuc"),   
          SIZE_BP_AGGR = 5,
          AGGREGATE_SAMPLES = FALSE,
          AGGREGATE_BIN = TRUE,
          MIN_FRAG_SIZE = 50,
          MAX_FRAG_SIZE =  250,
          AMPs = c("GAIN",  "LOSS")){
  
  # GET CNV REGIONS FROM PROGENETIX
  resCN <- getCNV_Regions(CLASS_CNV, FREQ_MANUAL = FREQ)
  
  res_Density <- parallel::mclapply(1:nrow(AllSample), function(i){
    
    sample <- AllSample[i,]
    print(sample$FASTQ_Name)
    
    load(sample$PathFragmentomics)
    
    df_ALT_GR <- resCN
    
    region_GR <- getRegionBinSample(sample$PathFeatures)
    
    ##new
    #names_region_GR <- paste(region_GR@seqnames, region_GR@ranges, sep = "-")
    
    print(table(df_ALT_GR$ALT))
    
    ALT_res <- lapply(AMPs, function(ALT){
      
      overlapGAIN <- findOverlaps(df_ALT_GR[df_ALT_GR$ALT==ALT,], region_GR, minoverlap = MIN_SIZE_ALT)
      
      ##new
      #indexOverlapTO <- unique(overlapGAIN@to)
      #resDff <- lapply(resFEATUREs[names_region_GR[indexOverlapTO]], function(x) x$region_weights)
      
      resDff <- lapply(resFEATUREs[unique(overlapGAIN@to)], function(x) x$region_weights) #old
    
      if(is.list(resDff)) resDff <- do.call(rbind, resDff)
      
      library(dplyr)
      resDff <- as.data.frame(resDff)
      resDff <- resDff %>%
        dplyr::group_by(Frag_len) %>%
        dplyr::summarise(Counts = sum(sum_GC_weight))
      resDff <- as.data.frame(resDff)
      resDff <- resDff[!is.na(resDff$Frag_len),]
      colnames(resDff) <- c("FragLen", "Counts")
      
      resDff
      
    })
    names(ALT_res) <- AMPs
    rm(resFEATUREs)
    
    ALT_res <- lapply(ALT_res, function(x){
      x <- x[x$FragLen>= MIN_FRAG_SIZE & x$FragLen<=MAX_FRAG_SIZE,]
    })
    
    mergeDF <- merge(ALT_res[[1]], ALT_res[[2]], by = "FragLen")
    
    bin <- seq(MIN_FRAG_SIZE, MAX_FRAG_SIZE,SIZE_BP_AGGR)
    
    SumBin_DF <- lapply(c(2,3), function(ind){
      
      Counts <- unlist(lapply(1:(length(bin)-1), function(i){
        counts <- sum(mergeDF[mergeDF$FragLen>=bin[i] & mergeDF$FragLen<bin[i+1],][,ind])
        counts
      }))
      
      Counts <- (Counts/sum(Counts))#*100
      Counts <- data.frame(Counts = Counts, Bin = paste(bin[-length(bin)], bin[-1]-1, sep = "-"))
      Counts$Class <- names(ALT_res)[ind-1]
      
      Counts
    })
    rm(ALT_res)
    
    mergeDF <- merge(SumBin_DF[[1]], SumBin_DF[[2]], by = "Bin")
    rm(SumBin_DF)
    
    mergeDF$Bin <- factor(mergeDF$Bin, levels = paste(bin[-length(bin)], bin[-1]-1, sep = "-"))
    mergeDF <- mergeDF[order(mergeDF$Bin),]
    
    #GLOBAL FEATURES
    P <- mergeDF$Counts.x
    Q <- mergeDF$Counts.y
    x <- rbind(P,Q)
    
    library(philentropy)
    # Calculate KL divergence
    mergeDF$KL_divergence <- KL(x, unit = 'log2')
    print(mergeDF$KL_divergence[1])
    
    mergeDF$Sum <- sum(abs(P-Q))
    mergeDF$Sd <- sd(P-Q)
    rm(P,Q)

    final <- data.frame(Sum = mergeDF$Sum[1], Sd = mergeDF$Sd[1] ,KL_divergence = (mergeDF$KL_divergence[1]), row.names = sample$FASTQ_Name)
    rm(mergeDF)
    final
    
  }, mc.cores = NUM_THREADS)
  MTX_density <- Reduce(rbind, res_Density)
  rm(res_Density)
  
  #Local Features
  resECDF_ALL <- parallel::mclapply(1:nrow(AllSample), function(i){
    
    MTX <- getMtxDiff_eCDF_Features_SINGLE_SAMP(NULL, resCN, pathFragm = AllSample$PathFragmentomics[i], pathFeat = AllSample$PathFeatures[i], features_sel = features_sel, MIN_SIZE_ALT = MIN_SIZE_ALT)
    rownames(MTX) <- AllSample$FASTQ_Name[i]
    
    MTX
  }, mc.cores = NUM_THREADS)
  MTX <- as.data.frame(Reduce(rbind, resECDF_ALL))
  rm(resECDF_ALL)
  
  print(head(MTX))
  print(head(MTX_density))
  
  list_features_mtx <- list(MTX, MTX_density)
  list_features_mtx <- list_features_mtx[unlist(lapply(list_features_mtx, function(x) !is.null(x)))]
  MTX <- Reduce(cbind, list_features_mtx) 
  
  if(MIN_MAX_NORM){
    MTX <- normalizeMatrixDataset(MTX, AllSample)
  }else{
    MTX
  } 
  
}



####### WORKFLOW #####

# remotes::install_github("progenetix/pgxRpi")
# library("pgxRpi")
# frequency <- pgxLoader(type="cnv_frequency", output ='pgxfreq',
#                        filters=c("NCIT:C3224"))
# pgxFreqplot(frequency)

# Define mapping of class to frequency and file path
CLASS_PARAMS_WGS <- list(
  Colon          = list(freq = 25,  file = "data/progenetix/Colorectal_Carcinoma_NCIT_C2955.tsv"),
  Lung           = list(freq = 25,  file = "data/progenetix/Lung_Carcinoma_NCIT_C4878.tsv"),
  Prostate       = list(freq = 25,  file = "data/progenetix/Prostate_Carcinoma_NCIT_C4863.tsv"),
  Urothelial     = list(freq = 22.5,file = "data/progenetix/Urothelial_Carcinoma_NCIT_C4030.tsv"),
  EwS            = list(freq = 10,  file = "data/progenetix/Ewing_Sarcoma_NCIT_C4817.tsv"),
  Mesotelioma    = list(freq = 15,  file = "data/progenetix/Malignant_Mesothelioma_NCIT_C4456.tsv"),
  Melanoma       = list(freq = 20,  file = "data/progenetix/Melanoma_NCIT_C3224.tsv"),
  Breast         = list(freq = 25,  file = "data/progenetix/Breast_Carcinoma_NCIT_C4872.tsv"),
  Gliobastoma    = list(freq = 35,  file = "data/progenetix/Gliobastoma_NCIT_C3058.tsv"),
  Astrocytoma    = list(freq = 25,  file = "data/progenetix/Astrocytoma_NCIT_C60781.tsv"),
  Oligodendroglioma = list(freq = 25, file = "data/progenetix/Oligodendroglioma_NCIT_C3288.tsv"),
  MPNST          = list(freq = 25,  file = "data/progenetix/MPNST_NCIT_C3798.tsv"),
  MM             = list(freq = 25,  file = "data/progenetix/Multiple_Myeloma_NCIT_C3242.tsv"),
  Pancreatic     = list(freq = 15,  file = "data/progenetix/Pancreatic_Adenocarcinoma_NCIT_C8294.tsv")
)

getFeatureBasedOnCNV <- function(AllSample, 
                                 CLASS_CNV, 
                                 NUM_THREADS = 30,
                                 features_sel = c("mean", "ratio_NucCor_Nuc" , "ratio_NucCorChrom_Nuc", "coverage","coverageNucCore", "coverageChrom","coverageNuc"))