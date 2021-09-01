#!/usr/bin/env Rscript
#
#args <- commandArgs(trailing = TRUE)
#sampleID <- args[1]

chromthripsis_detection <- function(SV_file, SCNV_file){
  intersect_standard_bedpe <- read.table(file = SV_file, header=TRUE)
  sum(intersect_standard_bedpe$SVTYPE == "INS")
  tmp <- intersect_standard_bedpe
  tmp$chrom1 <- as.character(tmp$chrom1)
  tmp$chrom2 <- as.character(tmp$chrom2)
  tmp$SVTYPE <- as.character(tmp$SVTYPE)
  #tmp[tmp$SVTYPE == "INS",]$strand1 <- "+"
  #tmp[tmp$SVTYPE == "INS",]$strand2 <- "-"
  
  if(sum(tmp$SVTYPE == "TRA_INV") !=0){
    tmp[tmp$SVTYPE == "TRA_INV",]$SVTYPE <- "TRA"
  }
  if(sum((tmp$strand1 == "+") & (tmp$strand2 == "+") & (tmp$SVTYPE == "INV")) != 0){
    tmp[(tmp$strand1 == "+") & (tmp$strand2 == "+") & (tmp$SVTYPE == "INV"),]$SVTYPE <- "h2hINV"
  }
  if(sum((tmp$strand1 == "-") & (tmp$strand2 == "-") & (tmp$SVTYPE == "INV")) != 0){
    tmp[(tmp$strand1 == "-") & (tmp$strand2 == "-") & (tmp$SVTYPE == "INV"),]$SVTYPE <- "t2tINV"
  }
  
  tmp <- tmp[tmp$chrom1!="chrY",]
  tmp <- tmp[tmp$chrom2!="chrY",]
  #tmp <- tmp[tmp$SVTYPE!="TRA",]
  tmp <- tmp[tmp$SVTYPE!="INS",]
  SV_data <- SVs(chrom1=as.character(substr(tmp$chrom1,4,nchar(tmp$chrom1))), pos1=as.numeric(tmp$pos1),
                 chrom2=as.character(substr(tmp$chrom2,4,nchar(tmp$chrom2))), pos2=as.numeric(tmp$pos2), 
                 SVtype=as.character(tmp$SVTYPE), 
                 strand1=as.character(tmp$strand1), 
                 strand2=as.character(tmp$strand2))
  
  SCNV <- read.table(SCNV_file, header = TRUE)
  SCNV$chromosome <- as.character(SCNV$chromosome)
  SCNV<- SCNV[SCNV$chromosome!="chrY",]
  
  d <- data.frame(chrom=as.character(substr(SCNV$chromosome,4, nchar(SCNV$chromosome))), 
                  start=SCNV$start,
                  end=SCNV$end,
                  total_cn=SCNV$cn)
  
  ################################################################
  # merging adjacent CNV with sample copy number value
  d$numeric_chrom <- as.character(d$chrom)
  d[d$numeric_chrom == "X",]$numeric_chrom <- "23"
  d$numeric_chrom <- as.numeric(d$numeric_chrom)
  d$total_cn <- as.numeric(d$total_cn)
  tmp <- group_by(d,Diff = (cumsum(c(1,diff(total_cn)) | c(1,diff(numeric_chrom)) !=0 )))
  dd <- summarise(tmp, numeric_chrom = c(numeric_chrom)[1], start = c(start)[1],end = c(end)[length(c(end))], total_cn = mean(total_cn))
  dd$chrom <- as.character(dd$numeric_chrom)
  dd[dd$chrom == "23",]$chrom <- "X"
  dd <- data.frame(chrom=dd$chrom, start=dd$start, end=dd$end, total_cn=dd$total_cn)
  d <- dd
  ################################################################
  CN_data <- CNVsegs(chrom=as.character(d$chrom), 
                     start=d$start,
                     end=d$end,
                     total_cn=d$total_cn)
  chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data)
  
  return(chromothripsis)
}

chromthripsis_high_conf <- function(chromothripsis_chromSummary, qvalue_threshold){
  p1 <- chromothripsis_chromSummary$pval_fragment_joins
  p2 <- chromothripsis_chromSummary$chr_breakpoint_enrichment
  p3 <- chromothripsis_chromSummary$pval_exp_cluster
  high_confidence_index1 <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 6 &
                                    chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states >=7 &
                                    p.adjust(p1,method = "fdr") > qvalue_threshold &
                                    ((p.adjust(p2,method = "fdr") <= qvalue_threshold & is.na(p3))|
                                       (is.na(p2) & p.adjust(p3,method = "fdr") <= qvalue_threshold)|
                                       (p.adjust(p2,method = "fdr") <= qvalue_threshold | p.adjust(p3,method = "fdr") <= qvalue_threshold)))
  
  high_confidence_index2 <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 3 &
                                    chromothripsis_chromSummary$number_TRA >=4 &
                                    chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states >=7 &
                                    p.adjust(p1,method = "fdr") > qvalue_threshold)
  
  low_confidence_index <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 6 &
                                  chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states %in% c(4,5,6) &
                                  p.adjust(p1,method = "fdr") > qvalue_threshold &
                                  ((p.adjust(p2,method = "fdr") <= qvalue_threshold & is.na(p3))|
                                     (is.na(p2) & p.adjust(p3,method = "fdr") <= qvalue_threshold)|
                                     (p.adjust(p2,method = "fdr") <= qvalue_threshold | p.adjust(p3,method = "fdr") <= qvalue_threshold)))
  return(list(high_confidence_index1, high_confidence_index2, low_confidence_index))
}

library(ShatterSeek)
### provide directory to store the shatterseek results; need MODIFY
directory <- "./" 
### provide directory of SV results (*“_Manta_GRIDSS_intersect_one_high_confidence.bedpe)
SV_directory <- "/scratch/gq19/tg2182/HYPER-DUP/MantaGRIDSS"
### provide directory to of CNV (*.-T.final.call.threshold.cns)
SCNV_file <- "/scratch/gq19/tg2182/HYPER-DUP/cnvkit_WJ"
### Run shatterseek
results <- c()
for(sampleID in c("KAL0104","KAL0106", "N0067", "SMU087", "UP2050", "UP2133")){ ### provide a set of sample ID; need MODIFY
    SV_file <- paste0(SV_directory, "/",sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe")  ### file name of SV 
    SCNV_file <- paste0(SCNV_directory,"/", sampleID,"-T.final.call.threshold",".cns") ### file name of copy number
    tmp <- chromthripsis_detection(SV_file, SCNV_file) ### Prepare the data and run ShatterSeek
    results <- rbind(results, data.frame(sampleID = sampleID, tmp@chromSummary))
}
### output of all chromothripsis regions for all samples
results$inter_other_chroms[results$inter_other_chroms == ""] <- NA
results$inter_other_chroms_coords_all[results$inter_other_chroms_coords_all == " "] <- NA
results$inter_other_chroms_coords_all <- gsub("\n",";",results$inter_other_chroms_coords_all)
write.table(results, paste0(directory,"df_chromothripsis_chromSummary.txt"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t") ### output of chromothripsis for all samples

### filtering to obtain high-confident chromothripsis regions
df_chromothripsis_chromSummary <- read.table(paste0(directory,"df_chromothripsis_chromSummary.txt", header = TRUE))
results <- c()
for(sampleID in c("KAL0104","KAL0106", "N0067", "SMU087", "UP2050", "UP2133")){ ### provide a set of sample ID; need MODIFY
    index <- chromthripsis_high_conf(df_chromothripsis_chromSummary, qvalue_threshold=0.2)
    results <- c()
    if(length(unique(c(index[[1]], index[[2]])))!=0){
      results <- rbind(results, data.frame(confidence = "high",
                                           df_chromothripsis_chromSummary[unique(c(index[[1]], index[[2]])),]))
    }
    if(length(unique(c(index[[3]])))!=0){
      results <- rbind(results, data.frame(confidence = "low",
                                           df_chromothripsis_chromSummary[unique(index[[3]]),]))
    }
}
### output of high-confident chromothripsis regions for all samples
write.table(results, paste0(directory,"df_confident_chromothripsis_chromSummary.txt"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t") 



