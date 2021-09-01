### Load packages required
library("ggplot2")
library(RColorBrewer)
library(ggpubr)

library(VariantAnnotation)
library(stringr)
library(stringi)
library(devtools)
library(StructuralVariantAnnotation)
library(vcfR)

### Load themes for plots
theme1 <-  theme(axis.text=element_text(size=12,face="bold"), 
                 axis.title=element_text(size=14,face="bold"), 
                 axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
                 plot.title = element_text(size=14),
                 legend.text = element_text(size=12,face="bold"),
                 #legend.title = element_text(size=12,face="bold"),
                 legend.title = element_blank(),
                 legend.position="top")

theme2 <-  theme(axis.text=element_text(size=12,face="bold"), 
                 axis.title=element_text(size=14,face="bold"), 
                 plot.title = element_text(size=14),
                 legend.text = element_text(size=12,face="bold"),
                 #legend.title = element_text(size=12,face="bold"),
                 legend.title = element_blank(),
                 legend.position="top")

theme6 <-  theme(axis.text=element_text(size=12,face="bold"), 
                 axis.title=element_text(size=14,face="bold"), 
                 plot.title = element_text(size=14),
                 legend.text = element_text(size=12,face="bold"),
                 legend.title = element_blank(),
                 legend.position="none")


ShinySoSV_prediction_figures <- function(Candidate_callers, df_prediction, performance, callset){
  if(callset == "individual"){
    df2 <- data.frame(sampleID = rep(df_prediction$sampleID, length(Candidate_callers)),
                      value = unlist(c(df_prediction[, grep("fit",colnames(df_prediction))])),
                      category = rep(Candidate_callers, each = nrow(df_prediction)))
    
  }
  if(callset %in% c("union","intersection")){
    combine_SV_SVcaller <- c()
    for (i in c(1:length(Candidate_callers))){
      combine_SV_SVcaller <- c(combine_SV_SVcaller, paste0(Candidate_callers[i], Candidate_callers[!(c(1:length(Candidate_callers)) %in% i)],c("Union","Intersect")[c("union", "intersection") %in% callset]))
    }
    df2 <- data.frame(sampleID = rep(df_prediction$sampleID, length(combine_SV_SVcaller)),
                      value = unlist(c(df_prediction[, grep("fit",colnames(df_prediction))])),
                      category = rep(combine_SV_SVcaller, each = nrow(df_prediction)))
  }
  
  assign(paste0("p",1), ggplot(data=df2, aes(x=sampleID, y = value, fill = category)) +
           geom_bar(stat="identity", position=position_dodge())+
           # scale_fill_manual(values=brewer.pal(12, "Paired")[c(2,4,6,8,10)],drop=F)+
           ylab(paste0("Predicted ", performance)) +
           scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
           theme1)
  pdf(file="./Shiny-SoSV_predictions_across_sample.pdf",width=15,height=5)
  print(p1)
  dev.off()
  
  p2 <- ggplot(df2, aes(x=category, y=value, group=category)) + 
    geom_boxplot(aes(fill=category)) + 
    labs(fill = "SV callers") + 
    ylab(paste0("Predicted ", performance)) +
    xlab("SV callers") +
    theme1
  
  pdf(file="./Shiny-SoSV_predictions_across_caller.pdf",width=15,height=5)
  print(p2)
  dev.off()
  
}

ShinySoSV_prediction <- function(Candidate_callers, newdata, performance, callset){
  model_name1 <- paste0(c("sen", "pre_off", "F1_score")[c("sensitivity", "precision", "F1_score") %in% performance],
                        c("", "_UnionIntersect", "_UnionIntersect")[c("individual", "union", "intersection") %in% callset])
  load(paste0("./Shiny-SoSV/data/","gam",model_name1,"_callers.RData"))
  
  model_name <- paste0(c("sen", "pre_off", "F1_score")[c("sensitivity", "precision", "F1_score") %in% performance])
  if(callset == "individual"){
    for(i in c(1: length(Candidate_callers))){
      assign(paste0("fit_",model_name,"_",Candidate_callers[i]), predict(eval(parse(text = paste0("gam", model_name, "_", Candidate_callers[i]))), newdata, type = "response",se.fit = T,unconditional = TRUE)$fit)
    }
    df_prediction <- data.frame(cbind(newdata, do.call("cbind", lapply(paste0("fit_",model_name,"_", Candidate_callers),function(s) eval(parse(text=s))))))
    colnames(df_prediction) <- c(colnames(newdata), paste0("fit_",performance,"_", Candidate_callers))  
  }
  
  if(callset %in% c("union","intersection")){
    combine_SV_SVcaller <- c()
    for (i in c(1:length(Candidate_callers))){
      combine_SV_SVcaller <- c(combine_SV_SVcaller, paste0(Candidate_callers[i], Candidate_callers[!(c(1:length(Candidate_callers)) %in% i)],c("Union","Intersect")[c("union", "intersection") %in% callset]))
    }
    for(i in c(1: length(combine_SV_SVcaller))){
      assign(paste0("fit_", model_name,"_", combine_SV_SVcaller[i]), predict(eval(parse(text = paste0("gam", model_name,"_", combine_SV_SVcaller[i]))), newdata, type = "response",se.fit = T,unconditional = TRUE)$fit)
    }
    df_prediction <- data.frame(cbind(newdata, do.call("cbind", lapply(paste0("fit_",model_name,"_", combine_SV_SVcaller),function(s) eval(parse(text=s))))))
    colnames(df_prediction) <- c(colnames(newdata), paste0("fit_", performance,"_", combine_SV_SVcaller)) 
  }
  write.csv(df_prediction,file = "./Shiny-SoSV_prediction.csv",row.names = FALSE)
  ShinySoSV_prediction_figures(Candidate_callers, df_prediction, performance, callset)
  return(df_prediction)
}


# ShinySoSV_prediction <- function(Candidate_callers, newdata){
#   for(i in c(1: length(Candidate_callers))){
#     assign(paste0("fit_","F1_score","_",Candidate_callers[i]), predict(eval(parse(text = paste0("gam","F1_score","_", Candidate_callers[i]))), newdata, type = "response",se.fit = T,unconditional = TRUE)$fit)
#   }
#   df_prediction <- data.frame(cbind(newdata, do.call("cbind", lapply(paste0("fit_","F1_score","_", Candidate_callers),function(s) eval(parse(text=s))))))
#   colnames(df_prediction) <- c(colnames(newdata), paste0("fit_","F1_score","_", Candidate_callers))                    
#   write.csv(df_prediction,file = "./Shiny-SoSV_prediction.csv",row.names = FALSE)
# }


# ShinySoSV_prediction_figures <- function(Candidate_callers, prediction_file){
#   df_prediction <- read.csv(prediction_file)
#   df2 <- data.frame(sampleID = rep(df_prediction$sampleID, length(Candidate_callers)),
#                     value = unlist(c(df_prediction[, grep("fit",colnames(df_prediction))])),
#                     category = rep(Candidate_callers, each = nrow(df_prediction)))
#   
#   assign(paste0("p",1), ggplot(data=df2, aes(x=sampleID, y = value, fill = category)) +
#                                                geom_bar(stat="identity", position=position_dodge())+
#                                                scale_fill_manual(values=brewer.pal(12, "Paired")[c(2,4,6,8,10)],drop=F)+
#                                                scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
#                                                theme1)
#   pdf(file="./Shiny-SoSV_predictions_across_sample.pdf",width=15,height=5)
#   print(p1)
#   dev.off()
#   
#   p2 <- ggplot(df2, aes(x=category, y=value, group=category)) + 
#                       geom_boxplot(aes(fill=category)) + 
#                       labs(fill = "SV callers") + 
#                       ylab("Predicted F1 score") +
#                       xlab("SV callers") +
#                       theme6
# 
#   pdf(file="./Shiny-SoSV_predictions_across_caller.pdf",width=5,height=5)
#   print(p2)
#   dev.off()
# }

Spectrum_SV_type <- function(input_SV_count, threshold_total, threshold_relative_freq){
  df <- data.frame(SVTYPE = rep(colnames(input_SV_count)[2:ncol(input_SV_count)], each = nrow(input_SV_count)),
                   sampleID = rep(input_SV_count$sampleID, (ncol(input_SV_count)-1)),
                   Count = as.vector(unlist(input_SV_count[,2:ncol(input_SV_count)])))
  
  assign(paste0("p",2), ggplot(data=df, aes(x=sampleID, y=Count, fill=SVTYPE)) +
           geom_bar(stat="identity")+
           scale_fill_manual(values = brewer.pal(12, "Set3"), drop=F)+
           theme1)
  pdf(file="./Spectrum_SV_across_sample.pdf", width=19, height=9)
  print(p2)
  dev.off()
  
  p1 <- ggplot(df, aes(x = SVTYPE, y=Count)) +
    geom_boxplot()+
    theme2
  
  N_total <- rowSums(input_SV_count[,2:ncol(input_SV_count)])
  if(missing(threshold_total)){
    threshold_total <- mean(N_total)
  }
  if(missing(threshold_relative_freq)){
    threshold_relative_freq <- 0.5
  }
  
  df2 <- cbind(df, 
               relative_freq = df$Count/(rep(N_total,(ncol(input_SV_count)-1))), 
               N_total = rep(N_total,(ncol(input_SV_count)-1)))
  p2 <- ggplot(df2, aes(x = SVTYPE, y=relative_freq)) +
    geom_boxplot()+
    ylab("Relative frequency")+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
    theme2
  
  figure1 <- ggarrange(p1, p2,
                       labels = c("A", "B"), heights = c(1,1),
                       ncol = 2, nrow = 1)
  
  pdf(file="./Spectrum_SV_across_SVTYPE.pdf",width=9,height=6)
  print(figure1)
  dev.off()
  
  hyper_SV <- df2[df2$relative_freq > threshold_relative_freq & df2$N_total > threshold_total,]
  return(hyper_SV)
}


SVTYPE_stat_generate <- function(bedpe){
  All_SVTYPE <- unique(bedpe$SVTYPE)
  for(i in c(1: length(All_SVTYPE))){
    assign(paste0("N_", All_SVTYPE[i]), sum(bedpe$SVTYPE == All_SVTYPE[i]))
  }
  STAT_bed <- data.frame(do.call("cbind", lapply(paste0("N_", All_SVTYPE),function(s) eval(parse(text=s)))))
  colnames(STAT_bed) <- All_SVTYPE
  return(STAT_bed)
}

Summary_SV_type <- function(All_sampleID, All_input_df_name){
  summary_results <- c()
  for(i in c(1:length(All_sampleID))){
    sampleID <- All_sampleID[i]
    SVTYPE_count <- SVTYPE_stat_generate(eval(parse(text = All_input_df_name[i])))
    summary_results <- rbind(summary_results, SVTYPE_count)
  }
  summary_results <- data.frame(sampleID =  All_sampleID, summary_results)
  return(summary_results)
}

Spectrum_SV_bin_generate <- function(All_sampleID, All_input_df_name, gene_file, bedtools_dir, directory){
  df_breakpoints <- Summary_SV_breakpoint(All_sampleID, All_input_df_name, gene_file, bedtools_dir, directory)
  df_bin_all <- Spectrum_SV_bin(df_breakpoints)
  df_bin_all_hotspots <- Spectrum_SV_bin_define_hotspot(df_bin_all) ##### add HOTSPOT information
  write.table(df_bin_all_hotspots,"./df_bin_all_hotspots.txt", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  return(df_bin_all_hotspots)
}

Summary_SV_breakpoint <- function(All_sampleID, All_input_df_name, gene_file, bedtools_dir, directory){
  ####Another method without keeping df_geneAnnotated files
  df_breakpoints <- c()
  for(i in c(1 : length(All_sampleID))){
    sampleID <- All_sampleID[i] 
    input_df_name <- All_input_df_name[i]
    
    bedpe_geneAnnotated <- SV_breakpoint_gene_annotation(sampleID, input_df_name, gene_file, bedtools_dir, directory)
    if(nrow(bedpe_geneAnnotated) != 0){
      df_breakpoints <- rbind(df_breakpoints,  data.frame(sampleID = sampleID, 
                                                          chrom = c(as.character(bedpe_geneAnnotated$chrom1), as.character(bedpe_geneAnnotated$chrom2)),
                                                          pos = c(bedpe_geneAnnotated$pos1, bedpe_geneAnnotated$pos2),
                                                          SVTYPE = c(as.character(bedpe_geneAnnotated$SVTYPE), as.character(bedpe_geneAnnotated$SVTYPE)),
                                                          overlap_gene = c(as.character(bedpe_geneAnnotated$pos1_overlap_gene), 
                                                                           as.character(bedpe_geneAnnotated$pos2_overlap_gene))))
    }
  }
  # ####One method keep df_geneAnnotated files
  # for(i in c(1 : length(All_sampleID))){
  #   sampleID <- All_sampleID[i]
  #   input_df_name <- All_input_df_name[i]
  #   bedpe_geneAnnotated <- SV_breakpoint_gene_annotation(All_sampleID, input_df_name, gene_file, bedtools_dir, directory)
  #   assign(paste0(sampleID, "_df_geneAnnotated"), data.frame(bedpe))
  # }
  # save(list = paste0(All_sampleID, "_df_geneAnnotated"), file = "./input_SV_bed_geneAnnotated.Rdata")
  # 
  # All_sampleID = paste0("sample_",c(1:100))
  # load("./input_SV_bed_geneAnnotated.Rdata")
  # All_input_df_name <- paste0(All_sampleID, "_df_geneAnnotated")
  # df_breakpoints <- c()
  # for(i in c(1 : length(All_sampleID))){
  #   sampleID <- All_sampleID[i] 
  #   input_df_name <- All_input_df_name[i]
  #   bedpe_geneAnnotated <- eval(parse(text = All_input_df_name[i]))
  #   if(nrow(bedpe_geneAnnotated) != 0){
  #     df_breakpoints <- rbind(df_breakpoints,  data.frame(sampleID = sampleID, 
  #                                                         chrom = c(as.character(bedpe_geneAnnotated$chrom1), as.character(bedpe_geneAnnotated$chrom2)),
  #                                                         pos = c(bedpe_geneAnnotated$pos1, bedpe_geneAnnotated$pos2),
  #                                                         SVTYPE = c(as.character(bedpe_geneAnnotated$SVTYPE), as.character(bedpe_geneAnnotated$SVTYPE)),
  #                                                         overlap_gene = c(as.character(bedpe_geneAnnotated$pos1_overlap_gene), 
  #                                                                          as.character(bedpe_geneAnnotated$pos2_overlap_gene))))
  #   }
  # }
  
  write.table(df_breakpoints,"df_breakpoints.txt", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  return(df_breakpoints)
}
SV_breakpoint_gene_annotation <- function(sampleID, input_df_name, gene_file, bedtools_dir, directory){
  gene_bed <- read.table(gene_file, header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
  bedpe <- eval(parse(text = input_df_name))
  if(nrow(bedpe) !=0){
    bedpe$chrom1 <- as.character(bedpe$chrom1)
    bedpe$chrom2 <- as.character(bedpe$chrom2)
    bedpe$ID <- as.character(bedpe$ID)
    bedpe$ID_mate <- as.character(bedpe$ID_mate)
    
    SV_bed <- data.frame(chrom = c(bedpe$chrom1, bedpe$chrom2),
                         start = c(bedpe$pos1-1, bedpe$pos2-1),
                         end = c(bedpe$pos1, bedpe$pos2),
                         ID = c(bedpe$ID, bedpe$ID_mate))
    
    write.table(SV_bed, paste0(directory,"SV_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    
    intersect_file <- paste0(directory,"SV_gene","_intersect.bed")
    system(paste(bedtools_dir,"intersect -a", paste0(directory,"SV_tmp.bed"),
                 "-b", gene_file,
                 "-wo >", intersect_file))
    
    intersect <- read.table(intersect_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    colnames(intersect) <- c("SV_chrom", "SV_start","SV_end","SV_ID",
                             colnames(gene_bed),"overlap")
    
    all_index <- match(SV_bed$ID, intersect$SV_ID)
    overlap_gene <- rep(NA, nrow(SV_bed))
    for(i in c(1:nrow(SV_bed))){
      if(!is.na(all_index[i])){
        overlap_gene[i] <- paste(unique(intersect[which(match(intersect$SV_ID, SV_bed$ID) == i),]$geneName), collapse= ",")
      }
    }
  }else{
    overlap_gene <- c()
  }
  
  bedpe_geneAnnotated <- cbind(bedpe, 
                               pos1_overlap_gene = overlap_gene[1:(length(overlap_gene)/2)],
                               pos2_overlap_gene = overlap_gene[((length(overlap_gene)/2)+1): length(overlap_gene)])
  return(bedpe_geneAnnotated)
}



Spectrum_SV_bin_gene <- function(df_bin_all_hotspots){
  df_hotspots <- df_bin_all[df_bin_all$is_hotspot,]
  hotspots <- unique(df_hotspots$chrom_bin_labels)
  percentage_bkpt_in_gene <- c()
  peak_gene <- c()
  percentage_bkpt_in_peak_gene <- c()
  hotspot_type <- c()
  for(i in c(1:length(hotspots))){
    tmp <- df_hotspots[df_hotspots$chrom_bin_labels == hotspots[i],]$overlap_gene
    percentage_bkpt_in_gene <- c(percentage_bkpt_in_gene, sum(!is.na(tmp))/length(tmp))
    if(sum(!is.na(tmp)) != 0){
      df_tmp <- data.frame(table(tmp))
      peak_gene <- c(peak_gene, paste0(df_tmp[,]$tmp, collapse = ";"))
      percentage_bkpt_in_peak_gene <- c(percentage_bkpt_in_peak_gene, paste0(df_tmp[,]$Freq/length(tmp), collapse = ";"))
      hotspot_type <- c(hotspot_type, c("breakpoint", "sample")[ifelse(sum(df_hotspots[df_hotspots$chrom_bin_labels == hotspots[i],]$is_hotspot_breakpoint)>0, 1,2)])
    }else{
      peak_gene <- c(peak_gene, NA)
      percentage_bkpt_in_peak_gene <- c(percentage_bkpt_in_peak_gene, NA)
      hotspot_type <- c(hotspot_type, NA)
    }
  }
  
  df_hotspots_info <- data.frame(hotspots, percentage_bkpt_in_gene, peak_gene, percentage_bkpt_in_peak_gene, hotspot_type)
  
  tmp <- df_hotspots_info[!is.na(df_hotspots_info$peak_gene),]
  hotspot_gene_df <- c()
  for(i in c(1:nrow(tmp))){
    hotspot_gene_df <- rbind(hotspot_gene_df, data.frame(gene = unlist(strsplit(tmp[i,]$peak_gene,";")),
                                                         percentage_bkpt = unlist(strsplit(tmp[i,]$percentage_bkpt_in_peak_gene,";")),
                                                         hotspot_type = unlist(strsplit(tmp[i,]$hotspot_type,";"))))
  }
  return(hotspot_gene_df)
}


Spectrum_SV_bin <- function(df_breakpoints){
  df_bin_all <- c()
  for(chrom in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")){
    df2 <- df_breakpoints[df_breakpoints$chrom %in% chrom,]
    # breaks = seq(0,250e6,1e5)
    breaks = seq(0, 250e6, 1e6)
    bin_labels <- cut(df2$pos, breaks=breaks, labels=c(1:(length(breaks)-1)))
    bin <- cut(df2$pos, breaks=breaks)
    df_bin <- cbind(bin_labels,bin, breaks = breaks[bin_labels],df2)
    
    bin_table1 <- data.frame(table(df_bin$bin_labels))
    
    df_sample <- df_bin[!(duplicated(paste0(df_bin$sampleID,"_", df_bin$bin_labels))),]
    bin_table2 <- data.frame(table(df_sample$bin_labels))
    
    df_bin <- cbind(df_bin, 
                    count_breakpoints = bin_table1[df_bin$bin_labels,]$Freq,
                    count_sample = bin_table2[df_bin$bin_labels,]$Freq)
    
    df_bin_all <- rbind(df_bin_all,df_bin)
  }
  
  df_bin_all <- data.frame(df_bin_all, chrom_bin_labels = paste0(df_bin_all$chrom,"_", df_bin_all$bin_labels))
  return(df_bin_all)
}



Spectrum_SV_bin_define_hotspot <- function(df_bin_all){
  ##### add HOTSPOT information
  df2 <- df_bin_all
  df3 <- df2
  df3 <- df3[!duplicated(df3$chrom_bin_labels),]
  SV_hotspots <- rbind(df3[df3$count_breakpoints>mean(df3$count_breakpoints)+3*sd(df3$count_breakpoints),],
                       df3[df3$count_sample>mean(df3$count_sample)+3*sd(df3$count_sample),])
  
  hotspots <- unique(SV_hotspots$chrom_bin_labels)
  df_bin_all_hotspot <- data.frame(df_bin_all,
                           is_hotspot_breakpoint = df_bin_all$chrom_bin_labels %in% df3[df3$count_breakpoints>mean(df3$count_breakpoints)+3*sd(df3$count_breakpoints),]$chrom_bin_labels,
                           is_hotspot_sample = df_bin_all$chrom_bin_labels %in% df3[df3$count_sample>mean(df3$count_sample)+3*sd(df3$count_sample),]$chrom_bin_labels,
                           is_hotspot = df_bin_all$chrom_bin_labels %in% hotspots)
  return(df_bin_all_hotspot)
}


Lumpy_somatic_PASS_filter <- function(vcf_file,Tumor_name, Normal_name, minSU,directory){
  setwd(directory)
  vcf <- readVcf(vcf_file, "hg38")
  somatic_vcf <- vcf[geno(vcf)$SU[,Normal_name]==0,]
  PASS_somatic_vcf <- somatic_vcf[geno(somatic_vcf)$SU[,Tumor_name]>=minSU,]
  return(PASS_somatic_vcf)
}

Lumpy_bedpe_generate <- function(vcf){
  info = info(vcf)
  gr <- rowRanges(vcf)
  
  idx <- !(sapply(info$SVLEN, length))
  info$SVLEN[idx] <- NA
  info$SVLEN <- unlist(info$SVLEN)
  
  idx <- !(sapply(info$END, length))
  info$END[idx] <- NA
  info$END <- unlist(info$END)
  
  idx <- !(sapply(info$MATEID, length))
  if (isEmpty(idx)){
    info$MATEID <- NA
  }else{
    info$MATEID[idx] <- NA
    info$MATEID <- unlist(info$MATEID)
  }
  
  if(is.null(info$EVENT)){
    info$EVENT <- NA
  }
  if(is.null(info$SVINSLEN)){
    info$SVINSLEN <- NA
  }
  #strand=info$STRANDS
  #strand_test <- strand@unlistData
  #strand@partitioning
  #idx <- !(sapply(info$STRANDS, length))
  #info$STRANDS[idx] <- NA
  #info$STRANDS <- unlist(info$STRANDS)
  ALT <- unlist(gr$ALT)
  tmp <- gsub("\\:.*",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom1 = as.character(seqnames(gr))
  chrom2 <- chrom1
  chrom2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  
  tmp <- gsub(".*\\:",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub("\\[.*",'', tmp)
  tmp2 <- gsub("\\].*",'', tmp1)
  pos2 <- info$END
  pos2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  
  bedpe <- data.frame(
    chrom1=seqnames(gr),
    pos1=start(gr),
    chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
    pos2=pos2,
    SVTYPE=info$SVTYPE,
    SVLEN=info$SVLEN,
    REF=gr$REF,
    ALT=ALT,
    #strand=info$STRANDS,
    ID=names(gr),
    score=gr$QUAL,
    FILTER=gr$FILTER,
    MATEID=info$MATEID,
    EVENT=info$EVENT, 
    stringsAsFactors = FALSE
  )
  
  return(bedpe)
}

Lumpy_standard_bedpe_generate <- function(bedpe, main_CHROM){ ###NOT FINISHED YET!!!! NEED to CHECK AGAIN!!!!
  if(main_CHROM==TRUE){
    bedpe <- bedpe[(bedpe$chrom1 %in% paste0("chr",c(seq(1,22,1),"X","Y")) ) & (bedpe$chrom2 %in% paste0("chr",c(seq(1,22,1),"X","Y"))),]
  }
  ### BNDs recorded twice so only keep BNDs recorded first in VCF file
  ### Another way to do is keeping both and calculate the mean position; or keep the one recorded as POS in VCF format file
  bedpe <- bedpe[is.na(match(bedpe$ID, bedpe$MATEID)) | (c(1:nrow(bedpe)) < match(bedpe$ID, bedpe$MATEID)),] ### either don't have mate (i.e. not BND) or present first as BND
  
  SV_index <- c(1:nrow(bedpe))
  ID_tmp <- bedpe$EVENT
  ID_tmp[is.na(bedpe$EVENT)] <- SV_index[is.na(bedpe$EVENT)]
  event_index <- match(ID_tmp,unique(ID_tmp))
  
  strand1 <- rep(NA,nrow(bedpe))
  strand2 <- rep(NA,nrow(bedpe))
  strand1[bedpe$SVTYPE == "DEL"] <- "+"
  strand2[bedpe$SVTYPE == "DEL"] <- "-"
  strand1[bedpe$SVTYPE == "DUP"] <- "-"
  strand2[bedpe$SVTYPE == "DUP"] <- "+"
  
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-"
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-"
  
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "+"
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "-"
  
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"
  
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "-"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "+"
  
 
  SVTYPE <- bedpe$SVTYPE
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA_INV"
  SVTYPE[(strand1 != strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  
  standard_bedpe <- data.frame(
    chrom1 = as.character(bedpe$chrom1),
    pos1 = as.integer(bedpe$pos1),
    chrom2 = as.character(bedpe$chrom2),
    pos2 = as.integer(bedpe$pos2),
    SVTYPE = SVTYPE,
    SVLEN = bedpe$SVLEN,
    strand1 = strand1,
    strand2 = strand2,
    ID = paste0("Lumpy_",SV_index,"_1_",event_index),
    ID_mate = paste0("Lumpy_",SV_index,"_2_",event_index),
    #ALT = bedpe$ALT,
    ID_caller = bedpe$ID,
    FILTER_caller = bedpe$FILTER,
    stringsAsFactors = FALSE)
  
  return(standard_bedpe)
}

manta_bedpe_generate <- function(vcf){
  info = info(vcf)
  gr <- rowRanges(vcf)
  
  idx <- !(sapply(info$SVLEN, length))
  info$SVLEN[idx] <- NA
  info$SVLEN <- unlist(info$SVLEN)
  
  idx <- !(sapply(info$MATEID, length))
  info$MATEID[idx] <- NA
  info$MATEID <- unlist(info$MATEID)
  
  idx <- !(sapply(info$SVINSLEN, length))
  info$SVINSLEN[idx] <- NA
  info$SVINSLEN <- unlist(info$SVINSLEN)
  
  idx <- !(sapply(info$SVINSSEQ, length))
  info$SVINSSEQ[idx] <- NA
  info$SVINSSEQ <- unlist(info$SVINSSEQ)
  
  idx <- !(sapply(info$LEFT_SVINSSEQ, length))
  info$LEFT_SVINSSEQ[idx] <- NA
  info$LEFT_SVINSSEQ <- unlist(info$LEFT_SVINSSEQ)
  
  idx <- !(sapply(info$RIGHT_SVINSSEQ, length))
  info$RIGHT_SVINSSEQ[idx] <- NA
  info$RIGHT_SVINSSEQ <- unlist(info$RIGHT_SVINSSEQ)
  
  idx <- !(sapply(info$EVENT, length))
  info$EVENT[idx] <- NA
  info$EVENT <- unlist(info$EVENT)
  
  ALT <- unlist(gr$ALT)
  tmp <- gsub("\\:.*",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom1 = as.character(seqnames(gr))
  chrom2 <- chrom1
  chrom2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  
  tmp <- gsub(".*\\:",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub("\\[.*",'', tmp)
  tmp2 <- gsub("\\].*",'', tmp1)
  pos2 <- info$END
  pos2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  
  bedpe <- data.frame(
    chrom1=chrom1,
    pos1=start(gr),
    chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
    pos2=pos2,
    SVTYPE=info$SVTYPE,
    SVLEN=info$SVLEN,
    REF=gr$REF,
    ALT=ALT,
    ID=names(gr),
    score=gr$QUAL,
    FILTER=gr$FILTER,
    MATEID=info$MATEID,
    EVENT=info$EVENT, ### for Manta version > 1.4.0, BNDs of simple inversion will have same EVENT ID
    SVINSLEN=info$SVINSLEN,
    SVINSSEQ=info$SVINSSEQ,
    LEFT_SVINSSEQ=info$LEFT_SVINSSEQ,
    RIGHT_SVINSSEQ=info$RIGHT_SVINSSEQ,
    stringsAsFactors = FALSE
  )
  return(bedpe)
}

Manta_standard_bedpe_generate <- function(bedpe, main_CHROM){
  if(main_CHROM==TRUE){
    bedpe <- bedpe[(bedpe$chrom1 %in% paste0("chr",c(seq(1,22,1),"X","Y")) ) & (bedpe$chrom2 %in% paste0("chr",c(seq(1,22,1),"X","Y"))),]
  }
  ### BNDs recorded twice so only keep BNDs recorded first in VCF file
  ### Another way to do is keeping both and calculate the mean position; or keep the one recorded as POS in VCF format file
  bedpe <- bedpe[is.na(match(bedpe$ID, bedpe$MATEID)) | (c(1:nrow(bedpe)) < match(bedpe$ID, bedpe$MATEID)),] ### either don't have mate (i.e. not BND) or present first as BND
  
  SV_index <- c(1:nrow(bedpe))
  ID_tmp <- bedpe$EVENT
  ID_tmp[is.na(bedpe$EVENT)] <- SV_index[is.na(bedpe$EVENT)]
  event_index <- match(ID_tmp,unique(ID_tmp))
  
  strand1 <- rep(NA,nrow(bedpe))
  strand2 <- rep(NA,nrow(bedpe))
  strand1[bedpe$SVTYPE == "DEL"] <- "+"
  strand2[bedpe$SVTYPE == "DEL"] <- "-"
  strand1[bedpe$SVTYPE == "DUP"] <- "-"
  strand2[bedpe$SVTYPE == "DUP"] <- "+"
  
  ###[p[t
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-" 
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-" 
  
  ###t[p[
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "+" 
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "-"
  
  ###t]p]
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"
  
  ###]p]t
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "-"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "+"
  
  # strand1[grepl('A]',substr(bedpe$ALT,1,2))|grepl('C]',substr(bedpe$ALT,1,2))|grepl('T]',substr(bedpe$ALT,1,2))|grepl('G]',substr(bedpe$ALT,1,2))] <- "+"
  # strand2[grepl('A]',substr(bedpe$ALT,1,2))|grepl('C]',substr(bedpe$ALT,1,2))|grepl('T]',substr(bedpe$ALT,1,2))|grepl('G]',substr(bedpe$ALT,1,2))] <- "+"
  # 
  # strand1[grepl('\\[A',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl('\\[C',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|
  #           grepl('\\[T',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl('\\[G',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))] <- "-"
  # strand2[grepl('\\[A',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl('\\[C',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|
  #           grepl('\\[T',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl('\\[G',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))] <- "-"
  # 
  # strand1[grepl('A\\[',substr(bedpe$ALT,1,2))|grepl('C\\[',substr(bedpe$ALT,1,2))|grepl('T\\[',substr(bedpe$ALT,1,2))|grepl('G\\[',substr(bedpe$ALT,1,2))] <- "+"
  # strand2[grepl('A\\[',substr(bedpe$ALT,1,2))|grepl('C\\[',substr(bedpe$ALT,1,2))|grepl('T\\[',substr(bedpe$ALT,1,2))|grepl('G\\[',substr(bedpe$ALT,1,2))] <- "-"
  # 
  # strand1[grepl(']A',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl(']C',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|
  #           grepl(']T',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl(']G',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))] <- "-"
  # strand2[grepl(']A',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl(']C',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|
  #           grepl(']T',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))|grepl(']G',substr(bedpe$ALT, nchar(bedpe$ALT)-2, nchar(bedpe$ALT)))] <- "+"
  # 
  SVTYPE <- bedpe$SVTYPE
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA_INV"
  SVTYPE[(strand1 != strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  
  standard_bedpe <- data.frame(
    chrom1 = as.character(bedpe$chrom1),
    pos1 = as.integer(bedpe$pos1),
    chrom2 = as.character(bedpe$chrom2),
    pos2 = as.integer(bedpe$pos2),
    SVTYPE = SVTYPE,
    SVLEN = bedpe$SVLEN,
    strand1 = strand1,
    strand2 = strand2,
    ID = paste0("Manta_",SV_index,"_1_",event_index),
    ID_mate = paste0("Manta_",SV_index,"_2_",event_index),
    #ALT = bedpe$ALT,
    ID_caller = bedpe$ID,
    FILTER_caller = bedpe$FILTER,
    stringsAsFactors = FALSE)
  
  return(standard_bedpe)
}

simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "BND", # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), "INV",
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}
GRIDSS_bedpe_generate <- function(vcf){
  gr <- breakpointRanges(vcf)
  info <- info(vcf)
  bedpe <- data.frame(
                chrom1=seqnames(gr),
                start1=start(gr) - 1,
                end1=end(gr),
                chrom2=seqnames(partner(gr)),
                start2=start(partner(gr)) - 1,
                end2=end(partner(gr)),
                name=names(gr),
                score=gr$QUAL,
                FILTER=gr$FILTER,
                strand1=strand(gr),
                strand2=strand(partner(gr)),
                SVTYPE=simpleEventType(gr),
                svlen=gr$svLen,
                inslen=gr$insLen,
                mate_name=names(partner(gr)),
                event=info$EVENT[rownames(info) %in% names(gr)],
                REF=gr$REF,
                ALT=gr$ALT, 
                stringsAsFactors = FALSE)
  return(bedpe)
}

GRIDSS_bedpe_generate_se <- function(vcf){
  segr <- breakendRanges(vcf)
  info <- info(vcf)
  if(length(segr) != 0){
    bedpe <- data.frame(
                        chrom1=seqnames(segr),
                        start1=start(segr) - 1,
                        end1=end(segr),
                        chrom2=seqnames(segr),
                        start2=start(segr),
                        end2=end(segr) + 1,
                        name=names(segr),
                        score=segr$QUAL,
                        FILTER=segr$FILTER,
                        strand1=strand(segr),
                        strand2=NA,
                        SVTYPE="BND",
                        svlen=segr$svLen,
                        inslen=segr$insLen,
                        mate_name=NA,
                        event=info$EVENT[rownames(info) %in% names(segr)],
                        REF=segr$REF,
                        ALT=segr$ALT, 
                        stringsAsFactors = FALSE)
  }else{
    bedpe <- data.frame(matrix(ncol = 18, nrow = 0))
    colnames(bedpe) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "FILTER", "strand1", "strand2",
                         "SVTYPE", "svlen", "inslen", "mate_name", "event", "REF", "ALT")
  }
  return(bedpe)
}

GRIDSS_standard_bedpe_generate <- function(bedpe, main_CHROM){
  if(main_CHROM==TRUE){
    bedpe <- bedpe[(bedpe$chrom1 %in% paste0("chr",c(seq(1,22,1),"X","Y")) ) & (bedpe$chrom2 %in% paste0("chr",c(seq(1,22,1),"X","Y"))),]
  }
  if(nrow(bedpe) != 0){
    if(!is.na(bedpe$mate_name)[1]){
      # For pair of breakends, SVs recorded twice so only keep SVs with ID with gridss.+o
      bedpe <- bedpe[str_detect(bedpe$name, "gridss.+o"),]
    }
    
    SVTYPE <- bedpe$SVTYPE
    SVTYPE[(bedpe$strand1 == bedpe$strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
    SVTYPE[(bedpe$strand1 == bedpe$strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA_INV"
    SVTYPE[(bedpe$strand1 != bedpe$strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
    
    SV_index <- c(1:nrow(bedpe))
    ID_tmp <- bedpe$event
    ID_tmp[is.na( bedpe$event)] <- SV_index[is.na( bedpe$event)]
    event_index <- match(ID_tmp,unique(ID_tmp))
    
    standard_bedpe <- data.frame( chrom1 = as.character(bedpe$chrom1),
                                  pos1 = as.integer((bedpe$start1 + bedpe$end1)/2),
                                  chrom2 = as.character(bedpe$chrom2),
                                  pos2 = as.integer((bedpe$start2 + bedpe$end2)/2),
                                  SVTYPE = SVTYPE,
                                  SVLEN = bedpe$svlen,
                                  strand1 = bedpe$strand1,
                                  strand2 = bedpe$strand2,
                                  ID = paste0("GRIDSS_",SV_index,"_1_",event_index),
                                  ID_mate = paste0("GRIDSS_",SV_index,"_2_",event_index),
                                  ID_caller = bedpe$name,
                                  FILTER_caller = bedpe$FILTER,
                                  stringsAsFactors = FALSE)
  }else{
    standard_bedpe <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(standard_bedpe) <- c("chrom1", "pos1", "chrom2", "pos2",  "SVTYPE", "SVLEN","strand1", "strand2",
                                  "ID", "ID_mate", "ID_caller", "FILTER_caller")
  }
  return(standard_bedpe)
}



SvABA_bedpe_generate <- function(vcf){
  SvABA_info = info(vcf)
  gr <- rowRanges(vcf)
  
  tmp <- gsub("\\:.*",'', unlist(rowRanges(vcf)$ALT))
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom2 <- tmp2
  
  SvABA_tmp <- data.frame(
    chrom1 = as.character(seqnames(rowRanges(vcf))),
    pos1 = start(rowRanges(vcf)),
    chrom2 = chrom2,
    pos2 = NA,
    ALT1 = as.character(unlist(rowRanges(vcf)$ALT)),
    ALT2 = NA,
    ID = names(gr),
    ID_mate = SvABA_info$MATEID,
    ID_event = gsub(":.*",'', names(gr)),
    ID_POS = gsub(".*:",'', names(gr)),
    SPAN = SvABA_info$SPAN,
    SVTYPE = SvABA_info$SVTYPE,
    QUAL=gr$QUAL,
    REF=gr$REF,
    FILTER=gr$FILTER,
    stringsAsFactors = FALSE
  ) 
  
  bedpe <- SvABA_tmp[SvABA_tmp$ID_POS==1,]
  bedpe[match(SvABA_tmp[SvABA_tmp$ID_POS==2,]$ID_event,bedpe$ID_event),]$pos2 <- SvABA_tmp[SvABA_tmp$ID_POS==2,]$pos1
  bedpe[match(SvABA_tmp[SvABA_tmp$ID_POS==2,]$ID_event,bedpe$ID_event),]$ALT2 <- SvABA_tmp[SvABA_tmp$ID_POS==2,]$ALT1
  
  bedpe$SVTYPE[(grepl('T\\[|A\\[|G\\[|C\\[|N\\[',bedpe$ALT1) & grepl('\\]T|\\]A|\\]G|\\]C|\\]N',bedpe$ALT2)) & bedpe$SPAN != -1] <- "DEL"
  bedpe$SVTYPE[(grepl('\\]T|\\]A|\\]G|\\]C|\\]N',bedpe$ALT1) & grepl('T\\[|A\\[|G\\[|C\\[|N\\[',bedpe$ALT2)) & bedpe$SPAN != -1] <- "DUP/INS"
  bedpe$SVTYPE[((grepl('T\\]|A\\]|G\\]|C\\]|N\\]',bedpe$ALT1) & grepl('T\\]|A\\]|G\\]|C\\]|N\\]',bedpe$ALT2))|
                     (grepl('\\[T|\\[A|\\[G|\\[C|\\[N',bedpe$ALT2) & grepl('\\[T|\\[A|\\[G|\\[C|\\[N',bedpe$ALT1))) &
                    bedpe$SPAN != -1] <- "INV"
  return(bedpe)
}

SvABA_standard_bedpe_generate <- function(bedpe, main_CHROM){
  if(main_CHROM==TRUE){
    bedpe <- bedpe[(bedpe$chrom1 %in% paste0("chr",c(seq(1,22,1),"X","Y")) ) & (bedpe$chrom2 %in% paste0("chr",c(seq(1,22,1),"X","Y"))),]
  }
  ### BNDs recorded twice so only keep BNDs recorded first in VCF file
  #bedpe <- bedpe[bedpe$ID_POS==1,]
  
  SV_index <- c(1:nrow(bedpe))
  ID_tmp <- bedpe$ID_event
  ID_tmp[is.na(bedpe$ID_event)] <- SV_index[is.na(bedpe$ID_event)]
  event_index <- match(ID_tmp,unique(ID_tmp))
  
  strand1 <- rep(NA,nrow(bedpe))
  strand2 <- rep(NA,nrow(bedpe))
  #strand1[bedpe$SVTYPE == "DEL"] <- "+"
  #strand2[bedpe$SVTYPE == "DEL"] <- "-"
  #strand1[bedpe$SVTYPE == "DUP"] <- "-"
  #strand2[bedpe$SVTYPE == "DUP"] <- "+"
  
  strand1[grepl('\\[', bedpe$ALT1) & substr(bedpe$ALT1,1,1) == "["] <- "-"
  strand2[grepl('\\[', bedpe$ALT1) & substr(bedpe$ALT1,1,1) == "["] <- "-"
  
  strand1[grepl('\\[', bedpe$ALT1) & substr(bedpe$ALT1,1,1) != "["] <- "+"
  strand2[grepl('\\[', bedpe$ALT1) & substr(bedpe$ALT1,1,1) != "["] <- "-"
  
  strand1[grepl(']', bedpe$ALT1) & substr(bedpe$ALT1,1,1) != "]"] <- "+"
  strand2[grepl(']', bedpe$ALT1) & substr(bedpe$ALT1,1,1) != "]"] <- "+"
  
  strand1[grepl(']', bedpe$ALT1) & substr(bedpe$ALT1,1,1) == "]"] <- "-"
  strand2[grepl(']', bedpe$ALT1) & substr(bedpe$ALT1,1,1) == "]"] <- "+"
  
  SVTYPE <- bedpe$SVTYPE
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA_INV"
  SVTYPE[(strand1 != strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  
  standard_bedpe <- data.frame(
    chrom1 = as.character(bedpe$chrom1),
    pos1 = as.integer(bedpe$pos1),
    chrom2 = as.character(bedpe$chrom2),
    pos2 = as.integer(bedpe$pos2),
    SVTYPE = SVTYPE,
    SVLEN = bedpe$SPAN,
    strand1 = strand1,
    strand2 = strand2,
    ID = paste0("SvABA_",SV_index,"_1_",event_index),
    ID_mate = paste0("SvABA_",SV_index,"_2_",event_index),
    ALT = bedpe$ALT1,
    ID_caller = bedpe$ID,
    FILTER_caller = bedpe$FILTER,
    stringsAsFactors = FALSE)
  
  return(standard_bedpe)
}







Standard_stat_generate <- function(bedpe){
  N_DEL <- sum(bedpe$SVTYPE == "DEL")
  N_DUP <- sum(bedpe$SVTYPE == "DUP")
  N_INV <- sum(bedpe$SVTYPE == "INV") ### for intra-chromosomal inversion
  N_INS <- sum(bedpe$SVTYPE == "INS")
  N_TRA <- sum(bedpe$SVTYPE == "TRA")  ### inter-chromosomal translocation
  N_TRA_INV <- sum(bedpe$SVTYPE == "TRA_INV") ### inter-chromosomal translocation with invertion signature at fusion junction
  
  N_sum <- N_DEL + N_DUP + N_INV + N_INS + N_TRA + N_TRA_INV
  
  STAT_bed <- c(N_sum,N_DEL,N_DUP,N_INV,N_INS,N_TRA,N_TRA_INV )
  names(STAT_bed) <- c("N_sum","N_DEL","N_DUP","N_INV","N_INS","N_TRA", "N_TRA_INV")
  return(STAT_bed)
}

Standard_SVLEN_stat_generate <- function(bedpe){
  SVLEN <- bedpe$SVLEN[!is.na(bedpe$SVLEN)]
  stat <- c(sum(abs(SVLEN) < 50),
            sum(abs(SVLEN) >= 50 & abs(SVLEN) < 100),
            sum(abs(SVLEN) >= 100 & abs(SVLEN) < 1000000),
            sum(abs(SVLEN) >= 1000000),
            sum(is.na(bedpe$SVLEN)))
  names(stat) <- c("<50","50-100","100-1000000",">=1000000","NO_SVLEN")
  return(stat)
}

Standard_bedtool_prepare_bkpt <- function(standard_bedpe, BND_diff){
  diff <- BND_diff/2
  standard_bed_tmp <- data.frame(
    chrom = c(standard_bedpe$chrom1, standard_bedpe$chrom2),
    start = c(standard_bedpe$pos1-diff, standard_bedpe$pos2-diff),
    end = c(standard_bedpe$pos1+diff, standard_bedpe$pos2+diff),
    SVTYPE = c(standard_bedpe$SVTYPE, standard_bedpe$SVTYPE),
    ID = c(as.character(standard_bedpe$ID), as.character(standard_bedpe$ID_mate)),
    ID_mate = c(as.character(standard_bedpe$ID_mate), as.character(standard_bedpe$ID)), stringsAsFactors = FALSE)
  if (sum(standard_bed_tmp$start<0)!=0){standard_bed_tmp[standard_bed_tmp$start<0,]$start <- 0}
  return(standard_bed_tmp)
}

TypePosfilter <- function(intersect_file, SVTYPE_ignore){
  intersect <- read.table(intersect_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(intersect) <- c("Caller1_CHROM", "Caller1_POS","Caller1_END","Caller1_SVTYPE","Caller1_ID","Caller1_ID_mate","Caller2",
                             "Caller2_CHROM", "Caller2_POS","Caller2_END","Caller2_SVTYPE","Caller2_ID","Caller2_ID_mate","overlap")
  
  
  if(SVTYPE_ignore){
    intersect_Typefilter <- intersect
  }else{
    intersect_Typefilter <- intersect[(intersect$Caller1_SVTYPE == intersect$Caller2_SVTYPE)| 
                                        grepl("BND",intersect$Caller1_SVTYPE) | grepl("BND",intersect$Caller2_SVTYPE),]
                                      
  }
  if(nrow(intersect_Typefilter)!=0){
    ### Check Caller_1 ID_mate also listed in Caller_1 ID
    tmp_index <- lapply(intersect_Typefilter$Caller1_ID_mate,function(x) which(x==intersect_Typefilter$Caller1_ID))
    ### Further check Caller_2 ID_mate also listed in Caller_2 ID, which matched with Caller_1 ID
    BND_ID_match <- vector(length=nrow(intersect_Typefilter))
    for (i in 1: nrow(intersect_Typefilter)){
      BND_ID_match[i] <- intersect_Typefilter$Caller2_ID_mate[i] %in% intersect_Typefilter$Caller2_ID[tmp_index[[i]]]
    }
    intersect_TypePosBNDfilter <- intersect_Typefilter[BND_ID_match,]
  }else{
    intersect_TypePosBNDfilter <- intersect_Typefilter
  }
  return(intersect_TypePosBNDfilter)
}

SVCaller_bed_newID_generate2 <- function(SVCaller_bed,SVCaller_name){
  if(nrow(SVCaller_bed) == 0){return( data.frame(matrix(ncol=0,nrow=0)))}
  colnames(SVCaller_bed) <- c("chrom1","pos1","chrom2","pos2","SVTYPE","SVCallerID","SVCallerID_mate")
  SV_index_tmp <- c(1:length(SVCaller_bed$SVCallerID))
  ID_tmp <- SVCaller_bed$SVCallerID
  SV_mate_index_tmp <- ifelse(is.na(match(SVCaller_bed$SVCallerID_mate,ID_tmp)),SV_index_tmp,match(SVCaller_bed$SVCallerID_mate,ID_tmp))
  SV_index <- ifelse(SV_index_tmp <= SV_mate_index_tmp,SV_index_tmp,SV_mate_index_tmp)
  mate1_index <- ifelse(duplicated(SV_index),"2","1")
  mate2_index <- ifelse(mate1_index=="1","2","1")
  event_index <- SV_index
  ID <- paste(SVCaller_name,"_",SV_index,"_",mate1_index,"_",event_index,sep="")
  ID_mate <- paste(SVCaller_name,"_",SV_index,"_",mate2_index,"_",event_index,sep="")
  SVCaller_bed_newID <- data.frame(cbind(SVCaller_bed[,c(1:5)],ID,ID_mate,SVCaller_bed[,c(6:7)]),stringsAsFactors = FALSE)
  return(SVCaller_bed_newID)
}

SVCaller_union_intersect_generate <- function(standard_bedpe_name,sampleID, SVCaller_name,BND_diff,bkpt_T_callers,directory,sub_directory, SVTYPE_ignore,bedtools_dir){
  ### Each bed, convert to bed_tmp and written to bed_tmp file
  for (i in 1:length(SVCaller_name)){
    assign(paste0(SVCaller_name[i],"_standard_bedpe"), eval(parse(text=standard_bedpe_name[i])))
    assign(paste0(SVCaller_name[i],"_bed_tmp"), Standard_bedtool_prepare_bkpt(eval(parse(text=paste0(SVCaller_name[i],"_standard_bedpe"))),BND_diff))
    setwd(directory)
    setwd(paste0("./",sub_directory,"/tmp"))
    write.table(eval(parse(text=paste0(SVCaller_name[i],"_bed_tmp"))), paste0(sampleID, "_", SVCaller_name[i],"_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  }
  
  ### Union set
  # generate bed and bed with new name, bed with new name bed_tmp and written to bed_tmp file
  SVCaller_name_all <- paste0(paste0(SVCaller_name,collapse = ""),"ALL")
  SVCaller_bed_all <- do.call("rbind", lapply(paste0(SVCaller_name,"_standard_bedpe"),function(s) eval(parse(text=s))))
  SVCaller_bed_all_newID <- SVCaller_bed_all
  #SVCaller_bed_all_newID <- SVCaller_bed_newID_generate2(SVCaller_bed_all,SVCaller_name_all)
  SVCaller_bed_all_newID_tmp <- Standard_bedtool_prepare_bkpt(SVCaller_bed_all_newID,BND_diff)
  write.table(SVCaller_bed_all_newID_tmp, paste0(sampleID, "_", SVCaller_name_all,"_tmp.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  
  ### Intersect set
  # bedtools intersect union set bed_tmp with all SV caller bed_tmp
  intersect_file <- paste0(sampleID, "_", "all_",paste0(SVCaller_name,collapse = "_"),"_intersect.bed")
  overlap_f <- (BND_diff - bkpt_T_callers)/BND_diff
  system(paste(bedtools_dir,"intersect -a", paste0(sampleID, "_", SVCaller_name_all,"_tmp.bed"),
               "-b", paste(paste0(sampleID, "_", SVCaller_name,"_tmp.bed"), collapse = " "),
               "-names",paste(SVCaller_name,collapse = " "), "-f",overlap_f, "-wo >", intersect_file))
  
  intersect_filter <- TypePosfilter(intersect_file, SVTYPE_ignore)
  ### remove ID in tmp bed, only the ID in original bed
  intersect_filter <- intersect_filter[intersect_filter$Caller1_ID %in% SVCaller_bed_all_newID$ID,]
  
  ############################################
  ### work for three sv caller overlapping ###
  data <- NULL
  for(i in 1: length(unique(intersect_filter$Caller1_ID))){
    tmp <- paste(intersect_filter[intersect_filter$Caller1_ID == unique(intersect_filter$Caller1_ID)[i],]$Caller2_ID)
    data <- rbind(data,
                  c(unique(intersect_filter$Caller1_ID)[i],
                    unlist(lapply(SVCaller_name,function(s)paste(tmp[grepl(s,tmp)],collapse = ",")))))
  }
  overlap_bed <- data.frame(data[,-1], stringsAsFactors = FALSE)
  colnames(overlap_bed) <- SVCaller_name
  overlap_bed <- unique(overlap_bed)
  overlap_index <- data.frame(overlap_bed !="")
  ############################################
  
  ### Union set (pick coordinate of the first SV caller)
  union1 <- paste(c(overlap_bed[overlap_index[,1],1]),collapse = ",")
  union2 <- paste(c(overlap_bed[overlap_index[,2] & (rowSums(overlap_index) != ncol(overlap_index)),2]),collapse = ",")
  SVCaller_bed_union1 <- eval(parse(text=paste0(SVCaller_name[1],"_standard_bedpe")))[as.character(eval(parse(text=paste0(SVCaller_name[1],"_standard_bedpe$ID")))) %in% strsplit(union1,split=",")[[1]],]
  SVCaller_bed_union2 <- eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe")))[as.character(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe$ID")))) %in% strsplit(union2,split=",")[[1]],]
  SVCaller_bed_union <- rbind(SVCaller_bed_union1,SVCaller_bed_union2)
  
  ### Add ID of overlapping caller in union set
  overlap_match_index <- lapply(SVCaller_bed_union$ID, function(x) which(grepl(x, overlap_bed[,1]))[1])
  overlap_match_index[!sapply(overlap_match_index,length)] <- NA
  ID_overlap_Caller <- overlap_bed[unlist(overlap_match_index),2]
  ID_overlap_Caller[ID_overlap_Caller == ""] <- NA
  SVCaller_bed_union <- cbind(SVCaller_bed_union, ID_overlap_Caller)

  # ### Add ID tag for overlapping caller in union set
  # overlap_caller_index <- lapply(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe")))$ID, function(x) which(grepl(x,ID_overlap_Caller))[1])
  # overlap_caller_index[!sapply(overlap_caller_index,length)] <- NA
  # tmp <- cbind(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe"))), overlap_index = unlist(overlap_caller_index))
  # ID_overlap_Caller2 <- c()
  # for(i in c(1:length(ID_overlap_Caller))){
  #   ID_overlap_Caller2 <- c(ID_overlap_Caller2, paste0(tmp[!is.na(tmp$overlap_index) & tmp$overlap_index == i,]$ID_caller,collapse = ","))
  #   
  # }
  # ID_overlap_Caller2[ID_overlap_Caller2 == ""] <- NA
  # ### Add FILTER tag for overlapping caller in union set
  # overlap_caller_index <- lapply(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe")))$ID, function(x) which(grepl(x,ID_overlap_Caller))[1])
  # overlap_caller_index[!sapply(overlap_caller_index,length)] <- NA
  # tmp <- cbind(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe"))), overlap_index = unlist(overlap_caller_index))
  # FILTER_overlap_Caller <- c()
  # for(i in c(1:length(ID_overlap_Caller))){
  #   FILTER_overlap_Caller <- c(FILTER_overlap_Caller, paste0(tmp[!is.na(tmp$overlap_index) & tmp$overlap_index == i,]$FILTER,collapse = ","))
  # }
  # FILTER_overlap_Caller[FILTER_overlap_Caller == ""] <- NA
  # 
  # SVCaller_bed_combine_all <- cbind(SVCaller_bed_union, FILTER_overlap_Caller, ID_overlap_Caller2)
  
  ID_overlap_Caller2 <- c()
  FILTER_overlap_Caller <- c()
  for(i in c(1:nrow(SVCaller_bed_union))){
    ID <- unlist(strsplit(as.character(SVCaller_bed_union$ID_overlap_Caller[i]),","))
    ID_overlap_Caller2 <- c(ID_overlap_Caller2, paste0(GRIDSS_standard_bedpe[GRIDSS_standard_bedpe$ID %in% ID,]$ID_caller, collapse = ","))
    FILTER_overlap_Caller <- c(FILTER_overlap_Caller, paste0(GRIDSS_standard_bedpe[GRIDSS_standard_bedpe$ID %in% ID,]$FILTER_caller, collapse = ","))
  }
  ID_overlap_Caller2[ID_overlap_Caller2 == ""] <- NA
  FILTER_overlap_Caller[FILTER_overlap_Caller == ""] <- NA
  SVCaller_bed_combine_all <- cbind(SVCaller_bed_union, FILTER_overlap_Caller, ID_overlap_Caller2)
  
  
  return (SVCaller_bed_combine_all)
  
  # ### Intersection set (pick coordinate of the first SV caller)
  # intersect <- paste(c(overlap_bed[rowSums(overlap_index) == ncol(overlap_index),1]),collapse = ",")
  # SVCaller_bed_intersect <- eval(parse(text=paste0(SVCaller_name[1],"_standard_bedpe")))[as.character(eval(parse(text=paste0(SVCaller_name[1],"_standard_bedpe$ID")))) %in% strsplit(intersect,split=",")[[1]],]
  # 
  # ### Add ID of overlapping caller
  # IntersectID <- SVCaller_bed_intersect$ID
  # overlap_match_index <- lapply(IntersectID, function(x) which(grepl(x, overlap_bed[,1]))[1])
  # ID_overlap_Caller <- overlap_bed[unlist(overlap_match_index),2]
  # ID_overlap_Caller[ID_overlap_Caller == ""] <- NA
  # 
  # ### Add FILTER tag for overlapping caller
  # overlap_caller_index <- lapply(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe")))$ID, function(x) which(grepl(x,ID_overlap_Caller))[1])
  # overlap_caller_index[!sapply(overlap_caller_index,length)] <- NA
  # tmp <- cbind(eval(parse(text=paste0(SVCaller_name[2],"_standard_bedpe"))), overlap_index = unlist(overlap_caller_index))
  # FILTER_overlap_Caller <- c()
  # for(i in c(1:length(ID_overlap_Caller))){
  #   FILTER_overlap_Caller <- c(FILTER_overlap_Caller, paste0(tmp[!is.na(tmp$overlap_index) & tmp$overlap_index == i,]$FILTER,collapse = ","))
  # }
  # FILTER_overlap_Caller[FILTER_overlap_Caller == ""] <- NA
  # ### Final intersect bedpe
  # SVCaller_bed_intersect <- cbind(SVCaller_bed_intersect, ID_overlap_Caller, FILTER_overlap_Caller)
  # 
  # List <- list(SVCaller_bed_union, SVCaller_bed_intersect)
  # return(List)
}

chromthripsis_detection <- function(SV_file, SCNV_file){
  intersect_standard_bedpe <- read.table(file = SV_file, header=TRUE)
  #sum(intersect_standard_bedpe$SVTYPE == "INS")
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
  
  #  SV_data <- data.frame(chrom1=as.character(substr(tmp$chrom1,4,nchar(tmp$chrom1))), pos1=as.numeric(tmp$pos1),
  #                        chrom2=as.character(substr(tmp$chrom2,4,nchar(tmp$chrom2))), pos2=as.numeric(tmp$pos2), 
  #                        SVtype=as.character(tmp$SVTYPE), 
  #                        strand1=as.character(tmp$strand1), 
  #                        strand2=as.character(tmp$strand2))
  # write.table(SV_data,"/Users/tingting/Desktop/work at home/HYPER-DUP/SV_data.txt", row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  #  
  
  SCNV <- read.table(SCNV_file, header = TRUE)
  SCNV$chromosome <- as.character(SCNV$chromosome)
  SCNV<- SCNV[SCNV$chromosome!="chrY",]
  
  d <- data.frame(chrom=as.character(substr(SCNV$chromosome,4, nchar(SCNV$chromosome))), 
                  start=SCNV$start,
                  end=SCNV$end,
                  total_cn=SCNV$cn)
 
  # #d <- d[c(1:10),]
  # #d <- rbind(d,c("2",1234,5678,1),c("3",1234,4567,1),c("X",1234,5678,1))
  # #d %>% arrange(total_cn) %>%
  # #  group_by(Diff = cumsum(c(1,diff(total_cn)) ==0)) %>%
  # #  summarise(chrom = paste0(chrom, collapse = "/"), start = paste0(start, collapse = "/"),end = paste0(end, collapse = "/"), Values = mean(total_cn)) %>%
  # #  ungroup() %>% select(-Diff)
  # 
  # #group_by(arrange(d,total_cn),Diff = cumsum(c(1,diff(total_cn)) ==0))
  # #group_by(d,Diff = cumsum(c(1,diff(total_cn)) !=0))
  # #group_by(d,Diff = (cumsum(c(1,diff(total_cn)) !=0 & cumsum(c(1,diff(chrom))))))
  
  # d$numeric_chrom <- as.character(d$chrom)
  # d[d$numeric_chrom == "X",]$numeric_chrom <- "23"
  # d$numeric_chrom <- as.numeric(d$numeric_chrom)
  # d$total_cn <- as.numeric(d$total_cn)
  # tmp <- group_by(d,Diff = (cumsum(c(1,diff(total_cn)) | c(1,diff(numeric_chrom)) !=0 )))
  # dd <- summarise(tmp, numeric_chrom = c(numeric_chrom)[1], start = c(start)[1],end = c(end)[length(c(end))], total_cn = mean(total_cn))
  # dd$chrom <- as.character(dd$numeric_chrom)
  # dd[dd$chrom == "23",]$chrom <- "X"
  # dd <- data.frame(chrom=dd$chrom, start=dd$start, end=dd$end, total_cn=dd$total_cn)
  # #write.table(dd,"/Users/tingting/Desktop/work at home/HYPER-DUP/dd.txt", row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  # d <- dd
  
  
  # #tmp1 <- which(d[c(1:nrow(d)-1),]$chrom == d[c(2:nrow(d)),]$chrom & d[c(1:nrow(d)-1),]$total_cn == d[c(2:nrow(d)),]$total_cn)
  # #tmp2 <- tmp1+1 %in% tmp1
  # #dd <- cbind(tmp1,tmp2)
  # #ddd <- data.frame(chrom=d[tmp1,]$chrom, start=)
  # 
  # #which(d[c(2:nrow(d)),]$chrom == d[c(1:nrow(d)-1),]$chrom & d[c(2:nrow(d)),]$total_cn != d[c(1:nrow(d)-1),]$total_cn)
  # 
  #write.table(d,"/Users/tingting/Desktop/work at home/HYPER-DUP/d.txt", row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  # dd <- d
  # dd$total_cn[dd$total_cn == 0] <- 150000
  # dd$total_cn[is.na(dd$total_cn)] <- 0
  # #library(GenomicRanges)
  # dd <- as(dd,"GRanges")
  # cov <- coverage(dd, weight = dd$total_cn)
  # dd1 <- as(cov,"GRanges")
  # dd1 <- as.data.frame(dd1)
  # dd1 <- dd1[dd1$score !=0,]
  # dd1 = dd1[,c(1,2,3,6)]
  # names(dd1) <- names(d)[1:4]
  # dd1$total_cn[dd1$total_cn == 150000] <- 0
  # d= dd1; rm(dd)
  # write.table(dd1,"/Users/tingting/Desktop/work at home/HYPER-DUP/dd1.txt", row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  # # SCNV <- d
  
  CN_data <- CNVsegs(chrom=as.character(d$chrom), 
                     start=d$start,
                     end=d$end,
                     total_cn=d$total_cn)
  chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data)
  
  return(chromothripsis)
}
chromthripsis_high_conf <- function(chromothripsis_chromSummary, qvalue_threshold){
  # pvalue_threshold <- 0.2
  # high_confidence_index1 <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 6 &
  #                                   chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states >=7 &
  #                                   chromothripsis_chromSummary$pval_fragment_joins > pvalue_threshold &
  #                                   ((chromothripsis_chromSummary$chr_breakpoint_enrichment < pvalue_threshold & is.na(chromothripsis_chromSummary$pval_exp_cluster))|
  #                                      (is.na(chromothripsis_chromSummary$chr_breakpoint_enrichment) & chromothripsis_chromSummary$pval_exp_cluster < pvalue_threshold)|
  #                                      (chromothripsis_chromSummary$chr_breakpoint_enrichment < pvalue_threshold | chromothripsis_chromSummary$pval_exp_cluster < pvalue_threshold)))
  # 
  # high_confidence_index2 <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 3 &
  #                                   chromothripsis_chromSummary$number_TRA >=4 &
  #                                   chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states >=7 &
  #                                   chromothripsis_chromSummary$pval_fragment_joins > pvalue_threshold)
  # 
  # low_confidence_index <- which(chromothripsis_chromSummary$number_DEL + chromothripsis_chromSummary$number_DUP + chromothripsis_chromSummary$number_h2hINV + chromothripsis_chromSummary$number_t2tINV >= 6 &
  #                                 chromothripsis_chromSummary$max_number_oscillating_CN_segments_2_states %in% c(4,5,6) &
  #                                 chromothripsis_chromSummary$pval_fragment_joins > pvalue_threshold &
  #                                 ((chromothripsis_chromSummary$chr_breakpoint_enrichment < pvalue_threshold & is.na(chromothripsis_chromSummary$pval_exp_cluster))|
  #                                    (is.na(chromothripsis_chromSummary$chr_breakpoint_enrichment) & chromothripsis_chromSummary$pval_exp_cluster < pvalue_threshold)|
  #                                    (chromothripsis_chromSummary$chr_breakpoint_enrichment < pvalue_threshold | chromothripsis_chromSummary$pval_exp_cluster < pvalue_threshold)))
  # 
  
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
