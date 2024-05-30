library(ggplot2)
library(patchwork)
library(tidyverse)

source('~/Desktop/orfeo_LTS/CNA/segmentation/claspy/utils.R')

results <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/results.tsv', sep = '\t')
real_bp <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/bps.tsv', sep = '\t')
true_bp <- get_true_bp_list(real_bp)


for (row in seq(1, nrow(results))){
  row <- 1
  tmp <- results %>% slice(row)
  tmp_res <- parse_res(tmp$res)
  
  true_res <- true_bp[[]]
  
  
  
}

