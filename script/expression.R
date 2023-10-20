# This script is used for making plots from Fantom5, based on hg19

library(tidyverse)
library(ggplot2)
library(data.table)

wd <- getwd()
setwd(wd)
sub <- "data"

# Functions
## get data frame of expression levels of target regions
get_df <- function(chr,TSS,strand,gene) {
  
  Location_subset <- Location %>%
    filter(chr == !!chr) %>%
    filter(strand == !!strand) %>%
    mutate(distance_to_TSS = !!TSS - start) %>%
    filter(abs(distance_to_TSS) < 1000) %>%
    arrange(abs(distance_to_TSS))
  Location_index <- Location_subset[1,1]
  target_gene_data <- table[table$X00Annotation %like% Location_index,]
  target_gene_data <- target_gene_data[,6:ncol(target_gene_data)]
  target_gene_data <- t(target_gene_data)
  df <- data.frame(sample = rownames(target_gene_data),
                   reads = target_gene_data[,1])
  
  df$gene <- gene
  return(df)
}

## get log transform data frame for target genes
get_logdf <- function(df) {
  df <- df %>%
    left_join(sample,by="sample") %>%
    mutate(log_reads = log10(reads+1))
  return(df)
}

input<- "hg19.cage_peak_phase1and2combined_counts_ann_decoded.osc.txt.gz.extract.tsv"
table <- read.delim(file.path(wd,sub,input))
names(table)[names(table) == "X00Annotation"] <- "Annotation"
table <- table %>%
  separate(Annotation,c("Annotation_tmp","strand"),sep=",",remove=FALSE) %>%
  separate(Annotation_tmp,c("chr","start","end"))

## export to bed and cooperate the hg38 bed information
hg19 <- table %>%
  select(Annotation) %>%
  separate(Annotation,c("chr","start","end"))
write.table(hg19,file.path(wd,sub,"hg19.bed"),sep="\t",
            row.names=F,col.names=F,quote=F)

## combine with hg38
## 26 peaks do not have hg38 coordinates
## UCSC coordinate lift +1 to the start site
hg38 <- read.delim(file.path(wd,sub,"hglft_genome_42464_d807c0.bed"),header=F)
hg38 <- hg38[,1:4]
colnames(hg38) <- c("hg38_chr","hg38_start","hg38_end","Annotation_nostrand")
hg38 <- hg38 %>%
  separate(Annotation_nostrand,c("chr","start","end")) %>%
  mutate(start = as.numeric(start) - 1) %>%
  mutate(index=paste(chr,start,sep="_")) %>%
  select(index,hg38_chr,hg38_start,hg38_end)

## Add hg38 information
table <- table %>%
  mutate(index=paste(chr,start,sep="_")) %>%
  left_join(hg38,by="index")

## re-organize
count <- cbind(table[,1:11],table[,90:92],table[,12:88])

# Make meta data
sample <- data.frame(sample=colnames(count[,15:ncol(count)]))
sample <- sample %>%
  mutate(sample2 = sample) %>%
  separate(sample2,c(NA,NA,NA,NA,NA,NA,"time","donor")) %>%
  separate(time,sep=4,c("hr","min")) %>%
  separate(hr,sep=2,c("hr",NA)) %>%
  separate(min,sep=2,c("min",NA)) %>%
  mutate(min = ifelse(min=="",0,min)) %>%
  mutate(min = as.numeric(min),hr = as.numeric(hr)) %>%
  mutate(total_time_min = hr*60 + min)

################################################
#### subset data and make plots (one gene) #####

## Set some variables here
## hg38 coordinate
## IL1Bp1 "chr2:112836770"
## IL1RN(+) TSS = "chr2:113127588"
## IL37(+) TSS = "chr2:113,670,548-113,676,458" 
## No CAGE expression data for IL37

## Subset data
## IL1B super TAD
start <- 112836770
end <- 113127588
count_subset <- count %>%
  filter(hg38_chr == "chr2") %>%
  filter(hg38_start >= !!start & hg38_start <= !!end)

## Get target count
# peak <- "p1@IL1B"
# target <- 112836770

# peak <- "p1@IL1RN"
# target <- 113117935

# peak <- "hg22100.1"
# target <-113124333

# peak <- "hg22093.1"
# target <- 113112256

peak <- "hg30274.1"
target <- 113105129

df <- count_subset %>%
  filter(hg38_start == !!target)

description <- df

## Combine data with meta data 
data <- data.frame(sample = colnames(df[,15:ncol(df)]),
                   count = t(df[,15:ncol(df)]))
data <- sample %>%
  left_join(data,by="sample") %>%
  group_by(total_time_min) %>%
  summarise(mean = mean(log10(count+1)),
            sd = sd(log10(count+1)))

## Make plots
p <- ggplot(data, aes(x=total_time_min, y=mean)) +
    geom_point() + geom_line() + theme_classic() +
  ylab("log10 RPM") +
  ggtitle(paste(peak,target,sep=" "))
p
fname <- file.path(wd,sub,paste("LPS_time_course_",peak,"_",target,".pdf"))
ggsave(fname)

################################################
################################################


################################################
######## Plot IL1B and IL1RN together ##########
# "p1@IL1B", "p1@IL1RN", "hg22100.1"
data <- data.frame()
for (target in c(112836770,113127588,113124333)) {
  df <- count_subset %>%
    filter(hg38_start == !!target)
  tmp <- data.frame(sample = colnames(df[,15:ncol(df)])
                    ,count = t(df[,15:ncol(df)]))
  tmp$short_description <- df$short_description
  data <- rbind(data,tmp)
}


# combine with meta data and make plots
data <- sample %>%
  left_join(data,by="sample") %>% 
  group_by(short_description,total_time_min) %>%
  summarise(mean = mean(log10(count+1)),
            sd = sd(log10(count+1)))
## Make plots
p <- ggplot(data, aes(x=total_time_min, y=mean, color=short_description)) +
  geom_point() + geom_line() +
  theme_classic() + ylab("log10 RPM")
p
fname <- file.path(wd,sub,paste("LPS_time_course_combined",".pdf"))
ggsave(fname)


## Make plots within 300 mins
time <- 400
data2 <- data %>%
  filter(total_time_min < !!time)

p2 <- ggplot(data2, aes(x=total_time_min, y=mean, color=short_description)) +
  geom_point() + geom_line() +
  theme_classic() + ylab("log10 RPM")
p2
fname <- file.path(wd,sub,paste("LPS_time_course_combined_",time,".pdf"))
ggsave(fname)

################################################
################################################