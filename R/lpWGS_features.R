
##### Download new CNV profile from Progenetix (TO TEST) ##### 

addNewCNV <- function(NCIT_CODE = "C3224"){
  
  OUTPUT_PATH <- paste0(PATH_INITIAL, "data/progenetix/NCIT_", NCIT_CODE, ".tsv")
  
  if(!file.exists(OUTPUT_PATH)){
    if(!"pgxRpi" %in% rownames(installed.packages())){
      remotes::install_github("progenetix/pgxRpi")
    }
    
    library("pgxRpi")
    frequency <- pgxLoader(type="cnv_frequency", output ='pgxfreq',
                           filters=c(paste0("NCIT:",NCIT_CODE)))
    #pgxFreqplot(frequency)
    
    #group_id	reference_name	start	end	gain_frequency	loss_frequency	no
    frequency_df <- as.data.frame(frequency)
    
    write.table(frequency_df, file = OUTPUT_PATH, sep = "\t")
  }
}

getDirFragm <- function(OUTPUT_DIR = "output/WGS/", FRAGM_DIR = "FRAGM_BIN/"){
  dirSave <- paste0(PATH_INITIAL, OUTPUT_DIR, FRAGM_DIR)
  if (!dir.exists(dirSave)) {
    dir.create(dirSave)
  }
  dirSave
}

getPathFragm <- function(sample, OUTPUT_DIR = "output/WGS/", FRAGM_DIR = "FRAGM_BIN/",SUFFIX_SAVE_FILE = "_res_frag_motif.RData"){
  dirSave <- getDirFragm(OUTPUT_DIR, FRAGM_DIR)
  path_output <- paste0(dirSave,sample,"_",as.integer(BIN_SIZE_WGS),SUFFIX_SAVE_FILE)
  path_output
}


##### extract fragm lenght and end-motif in 3Mb region (output _res_frag_motif.RData) #####
saveFragmBIN_fromBam <- function(sample, 
                                 bam,
                                 FASTA_FILE,
                                 PATH_SAMTOOLS,
                                 SUFFIX_BAM = "_recal.bam",
                                 PATH_OUTPUT_GC = "WGS_alignment/output_folder/GC_correction_output",
                                 SUFFIX_SAVE_FILE = "_res_frag_motif.RData",
                                 OUTPUT_DIR = "output/WGS/",
                                 FRAGM_DIR = "FRAGM_BIN/"){

  path_output <- getPathFragm(sample, OUTPUT_DIR, FRAGM_DIR, BIN_SIZE_WGS, SUFFIX_SAVE_FILE)
  
  BEDFILE <- paste0(PATH_INITIAL, "/data/genome_hg38_", as.integer(BIN_SIZE_WGS), ".bed")
  
  library(parallel)
  library(dplyr)
  
  GC_bias <- read.table(paste0(PATH_INITIAL, PATH_OUTPUT_GC, "/", sample, SUFFIX_BAM, "/", sample, SUFFIX_BAM, "_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz"), sep = "|")
  
  tmp_dir <- "/temp"
  awk_file_filter = paste0(PATH_INITIAL, "scripts/filter_wgs.awk")
  
  df_BED <- as.data.frame(read.csv(BEDFILE, sep = "\t"))
  regions <- apply(df_BED, 1, function(x) paste(gsub("chr","",x[1]), as.integer(x[2]), as.integer(x[3]) , sep = "-"))
  
  resFEATUREs <- mclapply(regions, function(region){
    
    print(region)
    
    region_data <-  strsplit(region, "-")[[1]]
    region_data[1] <- paste0("chr", region_data[1])
    position <- paste0(region_data[1],":",region_data[2],"-", region_data[3])
    
    watson <- paste0(PATH_SAMTOOLS, " view ", bam, " -f 99 ", position, " | awk -v MIN_MAPQ=", 
                     MAPQ, " -v MAX_FRAGMENT_LEN=", MAX_FRAG_LENGHT, 
                     " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]), 
                     " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", 
                     region_data[6], " -f ", awk_file_filter)
    
    crick <- paste0(PATH_SAMTOOLS, " view ", bam, " -f 163 ", 
                    position, " | awk -v MIN_MAPQ=", MAPQ, " -v MAX_FRAGMENT_LEN=", 
                    MAX_FRAG_LENGHT, " -v CHR=", region_data[1], 
                    " -v R_START=", as.numeric(region_data[2]), " -v R_END=", 
                    as.numeric(region_data[3]), " -v R_ID=", region_data[6], 
                    " -f ", awk_file_filter)
    
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
      watsonFRAG = watsonFRAG[watsonFRAG$V9>=MIN_FRAG_LENGHT,]
      watsonFRAG$V1 <- "watsonFRAG"
    }
    if (!is.null(crickFRAG)){
      crickFRAG = crickFRAG[crickFRAG$V9>=MIN_FRAG_LENGHT,]
      crickFRAG$V1 <- "crickFRAG"
    }
    
    if(is.null(crickFRAG) & is.null(watsonFRAG)){
      actualStrand <- NULL
    }else if(!is.null(crickFRAG) & !is.null(watsonFRAG)){
      actualStrand <- rbind(watsonFRAG,crickFRAG)
    }else if(!is.null(crickFRAG)){
      actualStrand <- crickFRAG
    }else if(!is.null(watsonFRAG)){
      actualStrand <- watsonFRAG
    }
    
    
    if(!is.null(actualStrand)){
      
      regions <- paste0(actualStrand$V6,":",actualStrand$V7,"-", actualStrand$V8-1)
      regions <- as.list(regions)
      d <- 1:length(regions)
      chunks <- split(d, ceiling(seq_along(d)/5000))
      
      resChunks <- lapply(chunks, function(chunk){
        
        func2 <-  paste(PATH_SAMTOOLS, "faidx ",FASTA_FILE, do.call(paste, regions[chunk]))
        MOTIFS <- system(func2, intern = TRUE)
        MOTIFS[grep("chr", MOTIFS)] <- "-"
        MOTIFS <- (MOTIFS[!grepl("chr", MOTIFS)])
        MOTIFS <- do.call(paste0, as.list(MOTIFS))
        MOTIFS <- substr(MOTIFS, 2, nchar(MOTIFS))
        MOTIFS <- strsplit(MOTIFS, split = "-")
        MOTIFS
      })
      
      resChunks <- unlist(resChunks)
      region_weights <-  lapply(resChunks, function(seq){
        fragm_len <- nchar(seq)
        GC_count <- (lengths(regmatches(seq, gregexpr("G", seq))))+(lengths(regmatches(seq, gregexpr("C", seq))))
        round(GC_bias[(fragm_len-MIN_FRAG_LENGHT)+1,GC_count+1],5)
      })
      
      GC_weight <- unlist(region_weights)
      actualStrand <- cbind(actualStrand, GC_weight)
      
      region_motifs <-  lapply(1:length(resChunks), function(num_seq){
        seq <- resChunks[[num_seq]]
        nameStrand <- actualStrand[num_seq,]$V1
        if(nameStrand=="watsonFRAG"){
          substr(seq, 1, 4)
        }else if(nameStrand=="crickFRAG"){
          substr(seq, nchar(seq)-3, nchar(seq))
        }
        
      })
      
      end_motif <- unlist(region_motifs)
      actualStrand <- cbind(actualStrand, end_motif)
      
      region_motifs <- as.data.frame(actualStrand %>% dplyr::group_by(end_motif) %>% dplyr::summarise(sum_weight = sum(GC_weight), count = n()))
      colnames(region_motifs) <- c("Motif", "sum_GC_weight", "Counts")
      
      region_weights <- as.data.frame(actualStrand %>% dplyr::group_by(V9) %>% dplyr::summarise(sum_weight = sum(GC_weight), count = n()))
      colnames(region_weights) <- c("Frag_len", "sum_GC_weight", "num_frag")
      
      list_res <- list(region_weights, region_motifs)
      names(list_res) <- c("region_weights","region_motifs")
      list_res
    }else{
      NULL
    }
    
  }, mc.cores = NUM_THREADS)
  
  names(resFEATUREs) <- regions
  
  save(resFEATUREs, file = path_output)

}

##### Compute metrics in 3Mb region ####

getDirMetrics <- function(OUTPUT_DIR = "output/WGS/", METRICS_DIR = "METRICS_BIN/"){
  dirSave <- paste0(PATH_INITIAL, OUTPUT_DIR, METRICS_DIR)
  if (!dir.exists(dirSave)) {
    dir.create(dirSave)
  }
  dirSave
}

getPathMetrics <- function(sample, OUTPUT_DIR = "output/WGS/", METRICS_DIR = "METRICS_BIN/", MOTIF = FALSE){
  dirSave <- getDirMetrics(OUTPUT_DIR, METRICS_DIR)
  ifelse(MOTIF, path_output <- paste0(dirSave,sample, "_motif_bin_",as.integer(BIN_SIZE_WGS),"_DF.RData"), path_output <- paste0(dirSave,sample, "_fragm_bin_",as.integer(BIN_SIZE_WGS),"_DF.RData"))
  path_output
}

saveMetricsBIN <- function(sample,
                           GC_CORR = TRUE,
                           OUTPUT_DIR = "output/WGS/",
                           FRAGM_DIR = "FRAGM_BIN/",
                           METRICS_DIR = "METRICS_BIN/"){

#dirRead <- paste0(PATH_INITIAL, OUTPUT_DIR, FRAGM_DIR)
#dirSave <- paste0(PATH_INITIAL, OUTPUT_DIR, METRICS_DIR)

# if (!dir.exists(dirSave)) {
#   dir.create(dirSave)
# }
#   
# dir.create(dirSave)
# setwd(dirSave)

MIN_NUCLEOSOME_CORE <- 140
MIN_CHROMATOSOME <- 160
MIN_NUCLEOSOME <- 171
MAX_NUCLEOSOME <- 240

coverage_nucleosome_core <- function(frag_lengths){
  sum(frag_lengths>=MIN_NUCLEOSOME_CORE  & frag_lengths<= MIN_CHROMATOSOME-1)
}  

coverage_chromatosome <- function(frag_lengths){
  sum(frag_lengths>=MIN_CHROMATOSOME  & frag_lengths<= MIN_NUCLEOSOME-1)
}

coverage_nucleosome <- function(frag_lengths){
  sum(frag_lengths>=MIN_NUCLEOSOME  & frag_lengths<= MAX_NUCLEOSOME)
}    

#path_fragm_data <- paste0(dirRead, sample, "_", as.integer(BIN_SIZE_WGS), "_res_frag_motif.RData")

path_fragm_data <- getPathFragm(sample, OUTPUT_DIR, FRAGM_DIR, BIN_SIZE_WGS)

path_output <- getPathMetrics(sample, OUTPUT_DIR, METRICS_DIR, BIN_SIZE_WGS)

#path_output <- paste0(dirSave,sample, "_fragm_bin_",as.integer(BIN_SIZE_WGS),"_DF.RData")

print(sample)

load(path_fragm_data)

resFEATUREs <- resFEATUREs[unlist(lapply(resFEATUREs, function(x) !is.null(x$region_weights)))]

df <- data.frame(row.names = names(resFEATUREs))

df$mean <- NA
df$coverage <- NA
df$coverageNucCore <- NA 
df$coverageChrom <- NA
df$coverageNuc <- NA

for(i in 1:length(resFEATUREs)){
  
  dff_all <- resFEATUREs[[i]]
  
  dff <- dff_all$region_weights
  if(nrow(dff)>0){
    
    if(GC_CORR){
      fragment_lengths <-  rep(dff$Frag_len,dff$sum_GC_weight)
    }else{
      fragment_lengths <-  rep(dff$Frag_len,dff$num_frag)
    }
    
    if(length(fragment_lengths)>0){
      
      df$mean[i] = mean(fragment_lengths)
      df$coverage[i] = sum(fragment_lengths)
      df$coverageNucCore[i] = coverage_nucleosome_core(fragment_lengths)
      df$coverageChrom[i] = coverage_chromatosome(fragment_lengths)
      df$coverageNuc[i] = coverage_nucleosome(fragment_lengths)
      
    }
  }
  
}

save(df, file = path_output)

path_output <- getPathMetrics(sample, OUTPUT_DIR, METRICS_DIR, BIN_SIZE_WGS, MOTIF = TRUE)

df <- data.frame(row.names = names(resFEATUREs))

df$CCCA <- NA
df$CCAG <- NA
df$CCTG <- NA
df$TAAA <- NA
df$AAAA <- NA
df$TTTT <- NA

for(i in 1:length(resFEATUREs)){
  
  dff_all <- resFEATUREs[[i]]
  
  dff <- dff_all$region_motifs
  
  resDff <- dff
  
  endMotif_all <- read.csv(paste0(PATH_INITIAL,"data/end-motif.csv"), sep = ";", header = FALSE)
  
  resDff <- resDff[resDff$Motif %in% endMotif_all,]
  
  if(sum(!endMotif_all %in% resDff$Motif)>0){
    notMot <- endMotif_all[!endMotif_all %in% resDff$Motif]
    resDffnotMot <- data.frame(Motif = as.integer(notMot), sum_GC_weight = 0 , Counts = 0)  
    resDff <- rbind(resDff, resDffnotMot)
  }
  
  if(!is.null(resDff)){
    library(dplyr)
    
    if(GC_CORR){
      resDff <- resDff %>%
        group_by(Motif) %>%
        summarise(Counts = sum(sum_GC_weight))
      resDff <- as.data.frame(resDff)
    }else{
      resDff <- resDff %>%
        group_by(Motif) %>%
        summarise(Counts = sum(Counts))
      resDff <- as.data.frame(resDff)
    }
    
    resDff$Density <- (resDff$Counts/sum(resDff$Counts))*100
  }
  
  df$CCCA[i] = resDff[resDff$Motif=="CCCA",]$Density
  df$CCAG[i] = resDff[resDff$Motif=="CCAG",]$Density
  df$CCTG[i] = resDff[resDff$Motif=="CCTG",]$Density
  df$TAAA[i] = resDff[resDff$Motif=="TAAA",]$Density
  df$AAAA[i] = resDff[resDff$Motif=="AAAA",]$Density
  df$TTTT[i] = resDff[resDff$Motif=="TTTT",]$Density
}

save(df, file = path_output)

}

#### pre-processing metrics ####

preProcessingFragmetomicsFeaturesAsList <- function(AllSampleRatioList, NA_VALUE = "zero", SCALE_FEAT = FALSE, MINMAX_FEAT = FALSE, MOTIF = FALSE){
  
  samplesName <- names(AllSampleRatioList)
  
  if(!MOTIF){
    AllSampleRatioList <- lapply(samplesName, function(sample){  
      
      x_gc <- AllSampleRatioList[[sample]]
      
      x <- c()
      
      x$mean <- x_gc$mean
      x$coverage <- x_gc$coverage
      x$coverageNucCore <- x_gc$coverageNucCore
      x$coverageChrom <- x_gc$coverageChrom
      x$coverageNuc <- x_gc$coverageNuc
      x$ratio_NucCor_Nuc <- x_gc$coverageNucCore/x_gc$coverageNuc
      x$ratio_NucCorChrom_Nuc <- (x_gc$coverageChrom+x_gc$coverageNucCore)/x_gc$coverageNuc
      
      if(any(is.infinite(x$ratio_NucCor_Nuc)) | any(is.infinite(x$ratio_NucCorChrom_Nuc)) | any(is.infinite(x$ratio)) ){
        x$ratio_NucCor_Nuc[is.infinite(x$ratio_NucCor_Nuc)] <- max(x$ratio_NucCor_Nuc[!is.infinite(x$ratio_NucCor_Nuc)], na.rm = TRUE)
        x$ratio_NucCorChrom_Nuc[is.infinite(x$ratio_NucCorChrom_Nuc)]  <- max(x$ratio_NucCorChrom_Nuc[!is.infinite(x$ratio_NucCorChrom_Nuc)], na.rm = TRUE)
      }
      
      x
    })
    names(AllSampleRatioList) <- samplesName
  }
  
  AllSampleRatioList <- lapply(samplesName, function(sample){  
    
    x <- AllSampleRatioList[[sample]]
    
    x <- lapply(x, function(y){
      
      if(NA_VALUE == "zero"){
        y[is.na(y)] <- 0
      }else{
        y[is.na(y)] <- mean(y, na.rm = TRUE)
      }
      
      if(SCALE_FEAT){
        y <- as.numeric(scale(y))
      }
      
      if(MINMAX_FEAT){
        y <- minMaxNorm(y, min = min(y), max = max(y))
      }
      
      y
      
    })
    
    
  })
  names(AllSampleRatioList) <- samplesName
  AllSampleRatioList
}

#### GET REGION OF BINs ####

getRegionBinSample <- function(pathMetrics){
  
  AllSampleDF <- get(load(pathMetrics))
  
  regionsOK <- rownames(AllSampleDF)
  
  regionsOK <- as.data.frame(matrix(unlist(strsplit(regionsOK, "-")), ncol=3, byrow=TRUE))
  
  colnames(regionsOK) <- c("Chr", "Start", "End")
  regionsOK$Chr <- as.integer(regionsOK$Chr)
  regionsOK$Start <- as.integer(regionsOK$Start)
  regionsOK$End <- as.integer(regionsOK$End)
  
  region_GR <- GenomicRanges::makeGRangesFromDataFrame(regionsOK, seqnames.field = "Chr", start.field = "Start", end.field = "End", keep.extra.columns = FALSE)
  region_GR
}

#### compute local feature ####

myECDFsingleValue <- function(x){
    x <- sort(x)
    n <- length(x)
    if (n < 1) 
      stop("'x' must have 1 or more non-missing values")
    vals <- unique(x)
    rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n, 
                      method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    out <- data.frame(x = knots(rval), cdf = rval(knots(rval)))
    out
}

#features_sel = c("ent","mean", "std","cv", "ratio_NucCor_Nuc" , "ratio_Chrom_Nuc","ratio_NucCorChrom_Nuc", "ratio", "coverage","coverageNucCore", "coverageChrom","coverageNuc")
getMtxDiff_eCDF_Features_SINGLE_SAMP <- function(CNV_regions, pathMetrics, features_sel = c("mean", "ratio_NucCor_Nuc" , "ratio_NucCorChrom_Nuc", "coverage","coverageNucCore", "coverageChrom","coverageNuc"), ECDF = TRUE, MIN_SIZE_ALT = 1500000){
  
  region_GR <- getRegionBinSample(pathMetrics)
  
  resSample <- list(get(load(pathMetrics)))
  names(resSample) <- pathMetrics
  resSample <- preProcessingFragmetomicsFeaturesAsList(resSample)[[1]]
  
  overlapGAIN <- findOverlaps(CNV_regions[CNV_regions$ALT=="GAIN",], region_GR, minoverlap = MIN_SIZE_ALT)
  overlapLOSS <- findOverlaps(CNV_regions[CNV_regions$ALT=="LOSS",], region_GR, minoverlap = MIN_SIZE_ALT)
  
  plotFEAT <- lapply(features_sel, function(feat){
    
    resGAINfeat <- list(resSample[[feat]][overlapGAIN@to])
    resGAINfeat <- as.data.frame(Reduce(rbind, resGAINfeat))
    
    resLOSSfeat <- list(resSample[[feat]][overlapLOSS@to])
    resLOSSfeat <- as.data.frame(Reduce(rbind, resLOSSfeat))
    
    if(grepl("coverage", feat)){
      normCov <- apply(resLOSSfeat, 2, median)
      resGAINfeat <- resGAINfeat/normCov
      resLOSSfeat <- resLOSSfeat/normCov
    }
    
    minHIST <- min(min(resGAINfeat),min(resLOSSfeat))
    maxHIST <- max(max(resGAINfeat),max(resLOSSfeat))
    
    if(!ECDF){
      
      seqBIN <- seq(minHIST,maxHIST, length.out = 100)
      
      res_histGAIN <- apply(resGAINfeat, 2, function(x){
        res_hist <- hist(as.numeric(x), breaks = seqBIN)
        res_hist$density
      })
      
      res_histLOSS <- apply(resLOSSfeat, 2, function(x){
        res_hist <- hist(as.numeric(x), breaks = seqBIN)
        res_hist$density
      })
      
    }else{
      
      seqBIN <- seq(minHIST,maxHIST, 0.01) 
      
      res_histGAIN <- apply(resGAINfeat, 2, function(x){
        x <- as.numeric(x)
        res_ecdf <- myECDFsingleValue(x)
        res_ecdf <- as.data.frame(approx(res_ecdf$x, res_ecdf$cdf, seqBIN, yleft = 0, yright = 1, ties = "ordered"))
        res_ecdf$y
      })
      
      res_histLOSS <- apply(resLOSSfeat, 2, function(x){
        x <- as.numeric(x)
        res_ecdf <- myECDFsingleValue(x)
        res_ecdf <- as.data.frame(approx(res_ecdf$x, res_ecdf$cdf, seqBIN, yleft = 0, yright = 1, ties = "ordered"))
        res_ecdf$y
      })
      
    }
    
    rm(seqBIN, resGAINfeat, resLOSSfeat)
    
    res_histDIFF <- res_histGAIN-res_histLOSS
    rm(res_histGAIN, res_histLOSS)
    res_histDIFF <- as.data.frame(res_histDIFF)
    
    df <- data.frame(Sum = sum(abs(res_histDIFF[,1])), Sd = sd(res_histDIFF[,1]))
    colnames(df) <- paste(colnames(df),feat, sep = "_")
    rm(res_histDIFF)
    df
    
  })
  
  Reduce(cbind, plotFEAT)
  
}

#### Extract CNV regions from Progenetix ####

getCNV_Regions <- function(CLASS, FREQ_MANUAL = NULL, FREQ_MANUAL_GAIN = NULL, FREQ_MANUAL_LOSS = NULL){
  
  library(GenomicRanges)
  
  if (!CLASS %in% names(CLASS_PARAMS_WGS)) {
    stop("Unknown CLASS: ", CLASS, ", specify CNV file!")
  }
  
  params <- CLASS_PARAMS_WGS[[CLASS]]
  
  FREQ <- params$freq
  CNV_FREQ <- read.table(paste0(PATH_INITIAL, "data/progenetix/NCIT_",params$NCIT,".tsv"), header = TRUE) #params$file 
  
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
  #CNV_FREQ$no <- NULL

  CNV_GR <- GenomicRanges::makeGRangesFromDataFrame(CNV_FREQ, seqnames.field = "reference_name", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

  CNV_GR <- unlist(GenomicRanges::reduce(GenomicRanges::split(CNV_GR, ~ALT)))
  CNV_GR$ALT <- names(CNV_GR)
  CNV_GR <- CNV_GR[CNV_GR@seqnames %in% 1:22,]
  
  CNV_GR
}

#### Extract WGS features based on CNV regions from Progenetix ####

getFeatureBasedOnCNV <- function(AllSample, 
          CLASS_CNV, 
          FREQ = NULL, 
          MIN_SIZE_ALT = 1500000, 
          features_sel = c("mean", "ratio_NucCor_Nuc" , "ratio_NucCorChrom_Nuc", "coverage","coverageNucCore", "coverageChrom","coverageNuc"),   
          SIZE_BP_AGGR = 5,
          AGGREGATE_SAMPLES = FALSE,
          AGGREGATE_BIN = TRUE,
          MIN_FRAG_LENGHT_SIZE = 50,
          MAX_FRAG_SIZE =  250,
          BIN_SIZE_WGS = 3000000,
          AMPs = c("GAIN",  "LOSS")){
  
  # GET CNV REGIONS FROM PROGENETIX
  CNV_regions <- getCNV_Regions(CLASS_CNV, FREQ_MANUAL = FREQ)
  
  #PathFragmentomics
  #PathMetrics  
  
  res_Density <- parallel::mclapply(AllSample, function(sample){
    

    #print(sample$FASTQ_Name)
    
 
    pathFragm <- getPathFragm(sample, BIN_SIZE_WGS = BIN_SIZE_WGS)
    load(pathFragm)
    
    pathMetrics <- getPathMetrics(sample, BIN_SIZE_WGS = BIN_SIZE_WGS)
    
    #load(sample$PathFragmentomics)
    region_GR <- getRegionBinSample(pathMetrics)
    #region_GR <- getRegionBinSample(sample$PathFeatures)
    
    print(table(CNV_regions$ALT))
    
    ALT_res <- lapply(AMPs, function(ALT){
      
      overlap_with_CNV <- findOverlaps(CNV_regions[CNV_regions$ALT==ALT,], region_GR, minoverlap = MIN_SIZE_ALT)
      
      resDff <- lapply(resFEATUREs[unique(overlap_with_CNV@to)], function(x) x$region_weights) 
    
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
      x <- x[x$FragLen>= MIN_FRAG_LENGHT_SIZE & x$FragLen<=MAX_FRAG_SIZE,]
    })
    
    mergeDF <- merge(ALT_res[[1]], ALT_res[[2]], by = "FragLen")
    
    bin <- seq(MIN_FRAG_LENGHT_SIZE, MAX_FRAG_SIZE,SIZE_BP_AGGR)
    
    SumBin_DF <- lapply(c(2,3), function(ind){
      
      Counts <- unlist(lapply(1:(length(bin)-1), function(i){
        counts <- sum(mergeDF[mergeDF$FragLen>=bin[i] & mergeDF$FragLen<bin[i+1],][,ind])
        counts
      }))
      
      Counts <- (Counts/sum(Counts))
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

    final <- data.frame(Sum = mergeDF$Sum[1], Sd = mergeDF$Sd[1] ,KL_divergence = (mergeDF$KL_divergence[1]), row.names = sample)
    rm(mergeDF)
    final
    
  }, mc.cores = NUM_THREADS)
  MTX_density <- Reduce(rbind, res_Density)
  rm(res_Density)
  
  #Local Features
  resECDF_ALL <- parallel::mclapply(AllSample, function(sample){
    
    pathMetrics <- getPathMetrics(sample, BIN_SIZE_WGS = BIN_SIZE_WGS)
    MTX <- getMtxDiff_eCDF_Features_SINGLE_SAMP(CNV_regions, pathMetrics = pathMetrics, features_sel = features_sel, MIN_SIZE_ALT = MIN_SIZE_ALT)
    rownames(MTX) <- sample
    
    #MTX <- getMtxDiff_eCDF_Features_SINGLE_SAMP(CNV_regions, pathMetrics = AllSample$PathFeatures[i], features_sel = features_sel, MIN_SIZE_ALT = MIN_SIZE_ALT)
    #rownames(MTX) <- AllSample$FASTQ_Name[i]
    
    MTX
  }, mc.cores = NUM_THREADS)
  MTX <- as.data.frame(Reduce(rbind, resECDF_ALL))
  rm(resECDF_ALL)
  
  print(head(MTX))
  print(head(MTX_density))
  
  list_features_mtx <- list(MTX, MTX_density)
  list_features_mtx <- list_features_mtx[unlist(lapply(list_features_mtx, function(x) !is.null(x)))]
  MTX <- Reduce(cbind, list_features_mtx) 
  
  MTX
}


