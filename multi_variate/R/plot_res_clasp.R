library(ggplot2)
library(tidyverse)

#results <- read.table('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/claspy/mean/results.tsv', sep = '\t', header = T)
results <- read.table('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/claspy/median/results_maf_median.tsv', sep = '\t', header = T)
data <- readRDS('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/data/smooth_data_meth.RDS') %>% mutate(order = seq(1, n()))

chr <- data %>% 
  group_by(CHROM) %>% 
  mutate(start = min(order),
         end = max(order)) %>% 
  select(CHROM, start, end) %>% 
  distinct(CHROM, start, end) %>% 
  mutate(half = (start + end)/2) 

chr_start_end <- chr %>% select(-half) %>% 
  pivot_longer(cols = c(start,end)) %>% pull(value)


true_bp <- data %>% 
  group_by(seg_id) %>% 
  mutate(start = min(order),
         end = max(order)) %>% 
  select(seg_id, start, end) %>% 
  distinct(seg_id, start, end) %>% 
  pivot_longer(cols = c(start,end)) %>% pull(value)

new_tmp <- c()
for (idx in seq(1, length(true_bp)-1)){
  if (as.integer(true_bp[idx]) == (as.integer(true_bp[idx+1]) -1)){
    x <- NULL
  } else {
    new_tmp <- c(new_tmp, as.integer(true_bp[idx]))
  }
}
new_tmp <- c(new_tmp, as.integer(true_bp[length(true_bp)]))
new_tmp <- new_tmp[2:(length(new_tmp)-1)]

source('~/Desktop/claspy/multi_variate/R/utils.R')
all_true <- c()
all_false <- c()

acceptance_wd = 10
for (row in seq(1, nrow(results))){
  tmp <- results %>% slice(row)
  
  pred_res <- parse_res(tmp$res)
  true_res <- new_tmp
  
  true <- 0
  for (real_bp in true_res){
    for (pred_bp in pred_res){
      tp <- pred_bp %in% seq(real_bp-acceptance_wd, real_bp+acceptance_wd)
      
      if (tp == TRUE){
        true <- true + 1
      }
      
    }
  }
  false <- length(pred_res) - true
  if (false < 0){
    false = 0
  }  
  
  true <- true / length(true_res)
  false <- false / length(pred_res)
  
  all_true <- c(all_true, true)
  all_false <- c(all_false, false)
}
results <- results %>% mutate(pred_corr = all_true,
                              pred_wrong = all_false) 

for (row in seq(1, nrow(results))){
  tmp <- results[row,]
  pred_res <- parse_res(tmp$res)
  if (length(pred_res) > 1){
    plt <- plot_claspy(data, tmp, type = 'median')  
  }
}

tmp <- results[6,]
type = 'median'


plot_claspy <- function(data, tmp, type = 'mean'){
  name <- paste0(tmp$mode, '_', tmp$thr, '_', tmp$ws)
  
  pred_res <- parse_res(tmp$res)
  dr <- paste0(type, '_dr')
  baf <- paste0(type, '_baf')
  maf <- paste0(type, '_meth')
  
  
  plt <- ggplot() +
    geom_point(aes(x = data$order, y = data[[dr]], color = data$CN), size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_color_manual(values = col_CN) +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylim(0,6) +
    ylab(dr) +
    xlab('pos') +
    
    ggplot() +
    geom_point(aes(x = data$order, y = data[[baf]], color = data$CN),size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_color_manual(values = col_CN) +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylim(0,1) +
    ylab(baf) +
    xlab('pos') +
    
    ggplot() +
    geom_point(aes(x = data$order, y = data$gt_AF, color = data$CN),size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_color_manual(values = col_CN) +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylab('VAF') +
    xlab('pos') +
    ylim(0,1) +
    
    ggplot() +
    geom_point(aes(x = data$order, y = data[[maf]], color = data$CN),size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_color_manual(values = col_CN) +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylim(0,1) +
    ylab(maf) +
    xlab('pos') +

    
    plot_layout(nrow = 4, guides = 'collect') + plot_annotation(title = name) & theme_bw()
    ggsave('~/Desktop/prova.png', dpi = 400, width = 12, height = 11, units = 'in', plot = plt)
  
    ggsave(paste0('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/claspy/median/claspy_',type, '_', name,'_maf.png'), dpi = 400, width = 12, height = 11, units = 'in', plot = plt)
}


col_CN <- list( ' ' = 'gray60',
                '0:0' = 'darkblue',
                '1:1' = ggplot2::alpha('forestgreen', .8),
                '1:0' = 'steelblue',
                '2:0' = 'plum4',
                '2:1' = ggplot2::alpha('orange', .8),
                '2:2' = 'firebrick',
                '3:0' = 'hotpink4',
                "3:1" = 'chocolate3',
                "3:2" = 'palegreen4', 
                "3:3" = 'plum4',
                '4:0' = 'salmon3',
                '4:1' = 'mediumpurple3',
                "4:2" = 'khaki2',
                '5:0' = 'palevioletred3',
                "5:2" = 'indianred1')
