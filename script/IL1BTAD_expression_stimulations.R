library(tidyverse)
library(ggplot2)

# dir ====
wd <- getwd()
input_dir <- file.path(wd,"data")
DE_dir <- file.path(input_dir,"DE_filtered")

# read files ==== 
sampleSheet <- read.csv(file.path(DE_dir,"sampleSheet.csv"))
sampleSheet <- sampleSheet[4:nrow(sampleSheet),2:ncol(sampleSheet)] # get rid of control group


IL1B_TAD <- read.csv(file.path(DE_dir,"IL1B_TAD.csv"))
transcript_index <- read.delim(file.path(DE_dir,"hg38.cage_peak_phase1and2combined_fair_counts_ann.osc.txt.gz.extract_transcriptSheet_with_hg38BED.tsv"))

count <- read.delim(file.path(DE_dir,"hg38.cage_peak_phase1and2combined_fair_counts_ann.osc.txt.gz.extract.tsv"))
colnames(count)[1] <- "Annotation"
count <- count[,c(1,8:ncol(count))]

count_annotated <- transcript_index %>%
  left_join(count,by="Annotation")

# IL1B TAD ====
count_annotated_IL1B <- count_annotated %>%
  filter(hg38_chr == "chr2") %>%
  filter(hg38_start > 112712617 & hg38_end < 113181216)

# subset data ====
target <- "IL1RN"
subset <- count_annotated_IL1B %>%
  filter(grepl(!!target,short_description))

description <- subset[nrow(subset),]

data_temp <- subset[nrow(subset),][,19:ncol(subset)] # get rid of control group
data = data.frame(sample = (rownames(t(data_temp))),count=t(data_temp))
colnames(data) <- c("sample","count")

data <- data %>%
  left_join(sampleSheet,by="sample") %>%
  group_by(treatment) %>%
  summarise(mean = mean(log10(count)),sd=sd(log10(count)))

p<- ggplot(data, aes(x=factor(data$treatment, levels = c("control","mock","glucan","BCG","Candida","Cryptococcus",
                                                         "streptococci","IFN","Salmonella","Trehalose","lipopolysaccharide")), y=mean, fill=treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=12),
        axis.text.y = element_text(size=12)) +
  labs(title=paste(description$short_description,description$hg38_start,sep=" "),x="treatment",y="log10 RPM")
p
fname <- paste(file.path(wd,paste0(paste(description$short_description,description$hg38_start,sep=" "),".png")))
ggsave(fname)
