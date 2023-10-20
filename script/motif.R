# This script is used for making plots from Fantom5, based on hg19

library(tidyverse)
library(ggplot2)
library(data.table)

wd <- getwd()
setwd(wd)
sub <- "data"

# set some variables
input <- "IL1Bp.csv"
strand = -1

table <- read.csv(file.path(wd,sub,input))

# Calculate motif frequency
frequency <- table %>%
  filter(strand == !!strand) %>%
  mutate(position = ifelse(abs(Start.Position) < 500, 500, 
                           ifelse(abs(Start.Position) < 1000, 1000, 2000))) %>%
  group_by(Transcription.Factor.Name, position) %>%
  summarize(n=n())