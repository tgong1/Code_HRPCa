#!/usr/bin/env Rscript
#
args <- commandArgs(trailing = TRUE)

#directory = "/scratch/gq19/tg2182/"
bedtools_dir <- "/scratch/gq19/tg2182/bedtools2/bin/bedtools"
source("/scratch/gq19/tg2182/HYPER-DUP/R_script/functions.R")

Sample_info_new <- read.csv(file = paste0(directory,"/HYPER-DUP/Sample_info_new.csv"))
Sample_info_new <- Sample_info_new[!(Sample_info_new$SampleID %in% c(#"10651-1042378","16599-1159140",
                                    "12333","N0031","N0055","KAL0068","KAL0095","TSH013")),]

directory <- "/scratch/gq19/tg2182/HYPER-DUP"
SVCaller_name <- c("Manta","GRIDSS")
sub_directory <- paste0(SVCaller_name,collapse = "")
BND_diff <- 2000
bkpt_T_callers <- 200
SVTYPE_ignore=FALSE

i <- as.numeric(args[1])
sampleID <- Sample_info_new$SampleID[i]
print(paste0("start R	script of ",sampleID))
vcf_file <- c(paste0("manta_",sampleID,".vcf"), paste0("gridss_",sampleID,".sv.vcf"))

### convert Manta unfiltered vcf to bedpe
assign(paste0(SVCaller_name[1], "_vcf"), readVcf(paste0(directory,"/",SVCaller_name[1],"/",SVCaller_name[1],"_somatic_VCF/raw_file/", vcf_file[1]), "hg38"))
assign(paste0(SVCaller_name[1], "_bedpe"), manta_bedpe_generate(eval(parse(text = paste0(SVCaller_name[1],"_vcf")))))
assign(paste0(SVCaller_name[1], "_standard_bedpe"), Manta_standard_bedpe_generate(eval(parse(text = paste0(SVCaller_name[1], "_bedpe"))), main_CHROM = TRUE))

### convert GRIDSS unfiltered vcf to bedpe
vcf <- readVcf(paste0(directory,"/",SVCaller_name[2],"/",SVCaller_name[2],"_VCF/", vcf_file[2]), "hg38")
somatic_vcf <- vcf[geno(vcf)$QUAL[,1] < 0.06 * rowRanges(vcf)$QUAL,]
assign(paste0(SVCaller_name[2], "_vcf"), somatic_vcf)
assign(paste0(SVCaller_name[2], "_bedpe"), GRIDSS_bedpe_generate(eval(parse(text = paste0(SVCaller_name[2],"_vcf")))))

### convert GRIDSS high confidence vcf to bedpe
GRIDSS_somatic_filter_vcf <- readVcf(paste0(directory,"/",SVCaller_name[2],"/",SVCaller_name[2],"_somatic_VCF/", "gridss_",sampleID,".sv.somatic.vcf"), "hg38")
GRIDSS_somatic_filter_bedpe <- GRIDSS_bedpe_generate(GRIDSS_somatic_filter_vcf)
GRIDSS_somatic_filter_standard_bedpe <- GRIDSS_standard_bedpe_generate(GRIDSS_somatic_filter_bedpe, main_CHROM = TRUE)

### add high-confidence SV into unfiltered bedpe if not in unfiltered bedpe
GRIDSS_bedpe <- rbind(GRIDSS_bedpe, GRIDSS_somatic_filter_bedpe[!(GRIDSS_somatic_filter_bedpe$name %in% GRIDSS_bedpe$name),])
assign(paste0(SVCaller_name[2], "_standard_bedpe"), GRIDSS_standard_bedpe_generate(eval(parse(text = paste0(SVCaller_name[2], "_bedpe"))), main_CHROM = TRUE))
standard_bedpe_name <- paste0(SVCaller_name, "_standard_bedpe")
SVCaller_bed_combine <- SVCaller_union_intersect_generate(standard_bedpe_name,sampleID,SVCaller_name,BND_diff,bkpt_T_callers,directory,sub_directory,SVTYPE_ignore,bedtools_dir)

### add column to indicate if this SV in GRIDSS high confidence somatic SV list
ID_overlap_Caller2 <- SVCaller_bed_combine$ID_overlap_Caller2
ID_overlap_Caller2 <- as.character(ID_overlap_Caller2)
SVCaller_bed_combine$ID_caller <- as.character(SVCaller_bed_combine$ID_caller)
ID_overlap_Caller2[grepl("GRIDSS",SVCaller_bed_combine$ID)] <- SVCaller_bed_combine[grepl("GRIDSS",SVCaller_bed_combine$ID),]$ID_caller

index <- lapply(ID_overlap_Caller2, function(x) which(stri_detect_fixed(x, GRIDSS_somatic_filter_standard_bedpe$ID_caller)))
High_Confidence_GRIDSS_index <- c()
for(i in c(1:length(index))){
  High_Confidence_GRIDSS_index <- c(High_Confidence_GRIDSS_index, paste0(index[[i]],collapse = ","))
}
High_Confidence_GRIDSS_index[High_Confidence_GRIDSS_index==""] <- NA
SVCaller_bedpe_combine_all <- cbind(SVCaller_bed_combine, High_Confidence_GRIDSS_index)

write.table(SVCaller_bedpe_combine_all, file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_combine_all_",bkpt_T_callers,"bp",".bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")

print(paste0("end R script of ",sampleID))


### union_intersection_callset.R
for(i in c(1:nrow(Sample_info_new))){
sampleID <- Sample_info_new$SampleID[i]
combine_all_standard_bedpe <- read.table(file = paste0(directory,"/",sub_directory,"/",sampleID, "_Manta_GRIDSS_combine_all.bedpe"), header=TRUE)

union_high_confidence_bedpe <- combine_all_standard_bedpe[(combine_all_standard_bedpe$FILTER_caller == "PASS" & grepl("Manta", combine_all_standard_bedpe$ID)) | 
                                           (!is.na(combine_all_standard_bedpe$High_Confidence_GRIDSS_index)),]
intersect_both_high_confidence_bedpe <- combine_all_standard_bedpe[combine_all_standard_bedpe$FILTER_caller == "PASS" & grepl("Manta", combine_all_standard_bedpe$ID) & 
                                                    (!is.na(combine_all_standard_bedpe$High_Confidence_GRIDSS_index)),]
#intersect_one_high_confidence_bedpe <- combine_all_standard_bedpe[!is.na(combine_all_standard_bedpe$ID_overlap_Caller), ] ### intersection set of all unfiltered SV
intersect_one_high_confidence_bedpe <- union_high_confidence_bedpe[!is.na(union_high_confidence_bedpe$ID_overlap_Caller), ]

write.table(union_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID,"_Manta_GRIDSS_union_both_high_confidence.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
write.table(intersect_both_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_both_high_confidence.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
write.table(intersect_one_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
}

### Count concordant SV calls 
union_PASS_SV_results <- c()
intersect_both_PASS_results <- c()
intersect_one_PASS_results <- c()
for(i in c(1:nrow(Sample_info_new))){
  sampleID <- Sample_info_new$SampleID[i]
  union_PASS <- read.table(file = paste0(directory,"/",sub_directory,"/", sampleID,"_Manta_GRIDSS_union.bedpe"), header=TRUE)
  intersect_both_PASS <- read.table(file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_both_high_confidence.bedpe"), header=TRUE)
  intersect_one_PASS <- read.table(file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe"), header=TRUE)
  
  union_PASS_SV_results <- rbind(union_PASS_SV_results, c(sampleID=as.character(sampleID), Standard_stat_generate(union_PASS)))
  intersect_both_PASS_results <- rbind(intersect_both_PASS_results, c(sampleID=as.character(sampleID), Standard_stat_generate(intersect_both_PASS)))
  intersect_one_PASS_results <- rbind(intersect_one_PASS_results, c(sampleID=as.character(sampleID), Standard_stat_generate(intersect_one_PASS)))
}
write.csv(union_PASS_SV_results,file = "./HRPCa_union_high_confidence_SV_results.csv",row.names = FALSE)
write.csv(intersect_both_PASS_results,file = "./HRPCa_intersect_both_high_confidence_results.csv",row.names = FALSE)
write.csv(intersect_one_PASS_results,file = "./HRPCa_intersect_one_high_confidence_results.csv",row.names = FALSE)





# #df <- c()
# for(i in c(1: nrow(Sample_info_new))){
#   sampleID <- Sample_info_new$SampleID[i]
#   combine_all_standard_bedpe <- read.table(file = paste0(directory,"/",sub_directory,"/",sampleID, "_Manta_GRIDSS_combine_all.bedpe"), header=TRUE)
#   
#   union_high_confidence_bedpe <- combine_all_standard_bedpe[(combine_all_standard_bedpe$FILTER_caller == "PASS" & grepl("Manta", combine_all_standard_bedpe$ID)) | 
#                                                               (!is.na(combine_all_standard_bedpe$High_Confidence_GRIDSS_index)),]
#   #intersect_both_high_confidence_bedpe <- combine_all_standard_bedpe[combine_all_standard_bedpe$FILTER_caller == "PASS" & grepl("Manta", combine_all_standard_bedpe$ID) & 
#   #                                                                     (!is.na(combine_all_standard_bedpe$High_Confidence_GRIDSS_index)),]
#   #intersect_one_high_confidence_bedpe <- combine_all_standard_bedpe[!is.na(combine_all_standard_bedpe$ID_overlap_Caller), ]
#   
#   intersect_one_high_confidence_bedpe <- union_high_confidence_bedpe[!is.na(union_high_confidence_bedpe$ID_overlap_Caller), ]
#   
#   #write.table(union_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID,"_Manta_GRIDSS_union.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
#   #write.table(intersect_both_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_both_high_confidence.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
#   write.table(intersect_one_high_confidence_bedpe, file = paste0(directory,"/",sub_directory,"/", sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
#   
#   #intersect_standard_bedpe <- read.table(file = paste0(directory,"/",sub_directory,"/",sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe"), header=TRUE)
#   #df <- rbind(df, c(sampleID = sampleID, count=nrow(intersect_standard_bedpe),
#   #                  count_both_low = sum(intersect_standard_bedpe$FILTER_caller != "PASS" & (is.na(intersect_standard_bedpe$High_Confidence_GRIDSS_index)))))
#   #intersect_standard_bedpe[intersect_standard_bedpe$FILTER_caller != "PASS" & (is.na(intersect_standard_bedpe$High_Confidence_GRIDSS_index)),]
# }
  


