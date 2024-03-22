# This script is transferred the ctcf hg19 map to ctcf hg38 version
## I only generate CTCF sites for IL1B TAD
## IL1B TAD in hg19
## chr2:113470194-113938793
library(tidyverse)

# path 
wd <- getwd()
output_dir <- file.path(wd,"data")
hg19_ctcf_file <- file.path(wd,"data","hg19.motifs.txt") 
hg38_ctcf_il1b_bed <- file.path(wd,"data","hg38_ctcf_il1b.bed")
# variable


# read ctcf raw data
## http://hicfiles.s3.amazonaws.com/internal/motifs/hg19.motifs.txt
hg19_ctcf <- read.delim(hg19_ctcf_file,sep="\t")

## clean
colnames(hg19_ctcf) <- c("pattern.name","chr","start",
                         "stop","strand","score","p.value","q.value",
                         "matched.sequence")
hg19_ctcf <- hg19_ctcf %>%
  arrange(chr) %>%
  mutate(chr = paste0("chr",chr)) %>%
  mutate(site.id = paste(chr,start,stop,strand,sep="_"))

# get CTCF sites from specific regions
## get CTCF sites from chr2
hg19_ctcf_chr2 <- hg19_ctcf %>%
  filter(chr == "chr2")

## get CTCF sites from IL1B TAD (a wider range)
## chr2:113,366,906-114,068,905
hg19_ctcf_il1b <- hg19_ctcf %>%
  filter(chr == "chr2") %>%
  filter(start > 113366906 & stop < 114068905)

## make bed files
df <- hg19_ctcf %>%
  select(chr,start,stop)
fname <- "hg19_ctcf.bed"
write.table(df,file.path(output_dir,fname),sep="\t",row.names=F,quote=F, col.names=F)

df <- hg19_ctcf_chr2 %>%
  select(chr,start,stop)
fname <- "hg19_ctcf_chr2.bed"
write.table(df,file.path(output_dir,fname),sep="\t",row.names=F, quote=F, col.names=F)

df <- hg19_ctcf_il1b %>%
  select(chr,start,stop)
fname <- "hg19_ctcf_il1b.bed"
write.table(df,file.path(output_dir,fname),sep="\t",row.names=F,quote=F, col.names=F)

# get hg38 ctcf bed files
hg38_ctcf_bed <- read.delim(hg38_ctcf_il1b_bed,sep="\t",header=F)

## clean
hg38_ctcf_bed <- hg38_ctcf_bed[,1:(ncol(hg38_ctcf_bed)-1)]
colnames(hg38_ctcf_bed) <- c("hg38_chr","hg38_start","hg38_stop","hg19_coordinates")

# Make the new ctcf file with hg38 coordinate
hg38_ctcf_il1b <- cbind(hg38_ctcf_bed,hg19_ctcf_il1b)

## write files
df <- hg38_ctcf_il1b
fname <- "hg38_ctcf_il1b.bed"
write.table(df,file.path(output_dir,fname),sep="\t",row.names=F,quote=F, col.names=F)


##
IL1RN_CTCF <- hg38_ctcf_il1b %>%
  filter(hg38_start > 113127598 & hg38_stop < 113138800) 
