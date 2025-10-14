##### extract fragm lenght and end-motif in 3Mb region (output _res_frag_motif.RData)####

args = commandArgs( trailingOnly = TRUE )
sample = args[1] 
bam = args[2]
dirSave = args[3]
BIN_SIZE = as.numeric(args[4])
NUM_THREADS = as.numeric(args[5])
FASTA_FILE = args[6]
PATH_SAMTOOLS = args[7]
MAPQ = as.numeric(args[8])
MAX_FRAG_LENGHT = as.numeric(args[9])
MIN_FRAG = 20

#SUFFIX_BAM <- "_recal"
SUFFIX_BAM <- gsub(".bam", "", args[10])
PATH_OUTPUT_GC <- "GC_correction_output"
SUFFIX_SAVE_FILE <- "_res_frag_motif.RData"

path_output <- paste0(dirSave,sample,"_",as.integer(BIN_SIZE),SUFFIX_SAVE_FILE)

dir.create(dirSave)

#source("../../../funcClassification.R")

PATH_INITIAL <- "/home3/adefalco/Fate-AI/"

BEDFILE <- paste0(PATH_INITIAL, "/acc_files/genome_", as.integer(BIN_SIZE), ".bed")
#BEDFILE <- paste0("./acc_files/genome_", as.integer(BIN_SIZE), ".bed")

library(parallel)
library(dplyr)

setwd(dirSave)
setwd("../")

GC_bias <- read.table(paste0(PATH_OUTPUT_GC, "/", sample,SUFFIX_BAM, "/", sample, SUFFIX_BAM, "_gc_weights_4simsMean.2IQRoutliersRemoved.2IgaussSmoothed.txt.gz"), sep = "|")

setwd("../")
print(getwd())

tmp_dir <- "/temp"
#awk_file_filter = "./Script/filter.awk"
#awk_file_stats = "./Script/stats.awk"

awk_file_filter = paste0(PATH_INITIAL, "/Script/filter.awk")
awk_file_stats = paste0(PATH_INITIAL, "/Script/stats.awk")

df_BED <- as.data.frame(read.csv(BEDFILE, sep = "\t"))
#load("~/ctDNAanalysis/newRatio/ICC02T_fragm_bin_3e+06_DF.RData")
#regions <- rownames(df)

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
    watsonFRAG = watsonFRAG[watsonFRAG$V9>=MIN_FRAG,]
    watsonFRAG$V1 <- "watsonFRAG"
  }
  if (!is.null(crickFRAG)){
    crickFRAG = crickFRAG[crickFRAG$V9>=MIN_FRAG,]
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
      round(GC_bias[(fragm_len-MIN_FRAG)+1,GC_count+1],5)
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