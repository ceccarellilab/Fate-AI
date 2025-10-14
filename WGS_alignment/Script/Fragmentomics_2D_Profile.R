##### Create fragmentomics profile ####

#source('~/ctDNAanalysis/funcClassification.R')

args = commandArgs( trailingOnly = TRUE )
sample_name = args[1] 
dirRead = args[2]
dirSave = args[3]
BIN_SIZE = as.numeric(args[4])
NUM_THREADS = as.numeric(args[5])
dirWorkflow = args[6]

GC_CORR = TRUE

dir.create(dirSave)
setwd(dirSave)

###### EXTRACT FETURES #####

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


# frag_ratio_dinucleosome <- function(frag_lengths){
#   # compute the ratio of short to long fragments
#   short_frags = sum(frag_lengths>=240  & frag_lengths<= 324)
#   long_frags = sum(frag_lengths>=325  & frag_lengths<= 400)
#   if(short_frags > 0 & long_frags > 0)
#     ratio = short_frags / long_frags
#   else
#     0
# }


frag_ratio <- function(frag_lengths){
  # compute the ratio of short to long fragments
  short_frags = sum(frag_lengths<= 120)
  long_frags = sum(frag_lengths>=140  & frag_lengths<= 250)
  if(short_frags > 0 & long_frags > 0)
    ratio = short_frags / long_frags
  else
    0
}

normalized_entropy <- function(frag_lengths, sample_name = "", i = "", bins = 50){
  #frag_lengths = frag_lengths[(frag_lengths>=100  & frag_lengths<= 200)]
  #png(paste(sample_name,i,"normalized_entropy.png",sep="_"), height=1080, width=1080, res=150) 
  #hist(frag_lengths, breaks = bins, plot = TRUE)
  #dev.off()
  
  histogram = hist(frag_lengths, breaks = bins, plot = FALSE)$counts
  pdf = histogram / sum(histogram)
  pdf = pdf[pdf!=0.0]
  moments = pdf * log(pdf) / log(length(pdf))
  -sum(moments)
}

path_fragm_data <- paste0(dirRead, sample_name, "_", as.integer(BIN_SIZE), "_res_frag_motif.RData")
path_output <- paste0(dirSave,sample_name, "_fragm_bin_",as.integer(BIN_SIZE),"_DF.RData")

print(sample_name)

load(path_fragm_data)

resFEATUREs <- resFEATUREs[unlist(lapply(resFEATUREs, function(x) !is.null(x$region_weights)))]

df <- data.frame(row.names = names(resFEATUREs))

df$ent <- NA
df$std <- NA
df$mean <- NA
df$cv <- NA
df$mad <- NA
df$coverage <- NA
df$coverageNucCore <- NA 
df$coverageChrom <- NA
df$coverageNuc <- NA

df$short_frags <- NA
df$long_frags <- NA
df$short_frags_din <- NA
df$long_frags_din <- NA

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
      
      df$ent[i] = normalized_entropy(fragment_lengths, sample_name, i)
      df$std[i] = sd(fragment_lengths)
      df$mean[i] = mean(fragment_lengths)
      df$cv[i] = df$std[i] / df$mean[i]
      df$mad[i] = median(abs(fragment_lengths - median(fragment_lengths)))
      
      df$coverage[i] = sum(fragment_lengths)
      df$coverageNucCore[i] = coverage_nucleosome_core(fragment_lengths)
      df$coverageChrom[i] = coverage_chromatosome(fragment_lengths)
      df$coverageNuc[i] = coverage_nucleosome(fragment_lengths)
      
      df$short_frags[i] = sum(fragment_lengths<= 120)
      df$long_frags[i] = sum(fragment_lengths>=140  & fragment_lengths<= 250)
      
      df$short_frags_din[i] =  sum(fragment_lengths>=240  & fragment_lengths<= 324)
      df$long_frags_din[i] = sum(fragment_lengths>=325  & fragment_lengths<= 400)
      
    }
  }
  
}

save(df, file = paste0(dirSave,sample_name, "_fragm_bin_",as.integer(BIN_SIZE),"_DF.RData"))


# if (file.exists(paste0(dirSave, sample_name, "_fragm_bin_",as.integer(BIN_SIZE),"_DF.RData"))) {
#   load(paste0(dirSave, sample_name, "_fragm_bin_",as.integer(BIN_SIZE),"_DF.RData"))
#   listFeatures <- as.list(df)
#   listFeatures
# }else{
#   print(paste0(sample_name,"_not_found"))
# }


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

  endMotif_all <- read.csv(paste0(dirWorkflow,"/acc_files/end-motif.csv"), sep = ";", header = FALSE)
  
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


save(df, file = paste0(dirSave,sample_name, "_motif_bin_",as.integer(BIN_SIZE),"_DF.RData"))


# if (file.exists(paste0(dirSave, sample_name, "_motif_bin_",as.integer(binSize),"_DF.RData"))) {
#   load(paste0(dirSave, sample_name, "_motif_bin_",as.integer(binSize),"_DF.RData"))
#   listFeatures <- as.list(df)
#   #names(listFeatures) <- c("ratio", "ratio_dinucl","ent", "std", "mean", "cv", "mad", "coverage")
#   listFeatures
# }else{
#   print(paste0(sample_name,"_not_found"))
# }


