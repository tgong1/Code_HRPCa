#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
sampleID <- args[1]


bedpe_to_ChainFinder_rr <- function(standard_bedpe, sampleID){
#  if(sum(standard_bedpe$SVTYPE == "INS")!=0){
#    standard_bedpe[standard_bedpe$SVTYPE == "INS",]$strand1 <- "+"
#    standard_bedpe[standard_bedpe$SVTYPE == "INS",]$strand2 <- "-"
#  }
#  standard_bedpe$chrom1 <- as.character(standard_bedpe$chrom1)
#  standard_bedpe$chrom2 <- as.character(standard_bedpe$chrom2)

  ChainFinder_rr <- data.frame(sample = sampleID, num = c(1:nrow(standard_bedpe)),
                               chr1 = substr(standard_bedpe$chrom1,4,nchar(standard_bedpe$chrom1)),
                               pos1 = standard_bedpe$pos1,
                               chr2 = substr(standard_bedpe$chrom2,4,nchar(standard_bedpe$chrom2)),
                               pos2 = standard_bedpe$pos2,
                               str1 = ifelse(standard_bedpe$strand1=="+",0,1),
                               str2 = ifelse(standard_bedpe$strand2=="+",0,1),
                               site1 = standard_bedpe$ID,
                               site2 = standard_bedpe$ID_mate)
  return(ChainFinder_rr)
}

library(rtracklayer)
directory <- "/scratch/gq19/tg2182/HYPER-DUP/ChainFinder/"
#wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
#gunzip hg38ToHg19.over.chain.gz
chain = import.chain(paste0(directory,"hg38ToHg19.over.chain"))

### Prepare input data for SV
intersect_standard_bedpe <- read.table(file = paste0(directory,"../MantaGRIDSS/",sampleID, "_Manta_GRIDSS_intersect_one_high_confidence.bedpe"), header=TRUE)
### exclude all DUPs
intersect_standard_bedpe <- intersect_standard_bedpe[intersect_standard_bedpe$SVTYPE != "DUP",]

SV_pos1 <- data.frame(chrom = intersect_standard_bedpe$chrom1,
                      start = intersect_standard_bedpe$pos1-1,
                      end = intersect_standard_bedpe$pos1,
                      strand = intersect_standard_bedpe$strand1,
                      name = intersect_standard_bedpe$ID)

SV_pos2 <- data.frame(chrom = intersect_standard_bedpe$chrom2,
                 start = intersect_standard_bedpe$pos2-1,
                 end = intersect_standard_bedpe$pos2,
                 strand = intersect_standard_bedpe$strand2,
                 name = intersect_standard_bedpe$ID)

Hg38toHg19 <- function(df38, chain){
  gr38 <- makeGRangesFromDataFrame(df38, keep.extra.columns=TRUE)
  seqlevelsStyle(gr38) = "UCSC"  # necessary
  gr19 = liftOver(gr38, chain)
  gr19 = unlist(gr19)
  genome(gr19) = "hg19"
  df19 <- as.data.frame(gr19)
  return(df19)
}
SV_pos1_hg19 <- Hg38toHg19(SV_pos1, chain)
SV_pos2_hg19 <- Hg38toHg19(SV_pos2, chain)
SV_pos1_hg19$seqnames <- as.character(SV_pos1_hg19$seqnames)
SV_pos2_hg19$seqnames <- as.character(SV_pos2_hg19$seqnames)

ifelse(nrow(SV_pos1_hg19) <= nrow(SV_pos2_hg19), tmp1 <- SV_pos1_hg19, tmp1 <- SV_pos2_hg19)
ifelse(nrow(SV_pos1_hg19) <= nrow(SV_pos2_hg19), tmp2 <- SV_pos2_hg19, tmp2 <- SV_pos1_hg19)

standard_bedpe <- intersect_standard_bedpe[intersect_standard_bedpe$ID %in% tmp1$name,]
standard_bedpe$chrom1 <- tmp1$seqnames
standard_bedpe$pos1 <- tmp1$end
standard_bedpe$chrom2 <- tmp2[match(tmp1$name,tmp2$name),]$seqnames
standard_bedpe$pos2 <- tmp2[match(tmp1$name,tmp2$name),]$end

ChainFinder_rr <- bedpe_to_ChainFinder_rr(standard_bedpe, sampleID)

write.table(ChainFinder_rr,file = paste0(directory,sampleID,"_rr_hg19.txt"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")

setwd(paste0(directory, "/../cnvkit_WJ"))
SCNV <- read.table(paste0(sampleID,"-T.final",".cns"), header = TRUE)
SCNV_hg38 <- cbind(SCNV[,c(1:3)],name = c(1:nrow(SCNV)))
SCNV_hg19 <- Hg38toHg19(SCNV_hg38, chain)
SCNV_hg19$seqnames <- as.character(SCNV_hg19$seqnames)

ChainFinder_cn <- data.frame(sample = sampleID, 
                             chr = substr(SCNV_hg19$seqnames,4, nchar(SCNV_hg19$seqnames)),
                             start = SCNV_hg19$start,
                             end = SCNV_hg19$end,
                             num_probes = SCNV[SCNV_hg19$name,]$probes,
                             segment_mean = SCNV[SCNV_hg19$name,]$log2)
write.table(ChainFinder_cn,file = paste0(directory,sampleID,"_cn_hg19.txt"), row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")

