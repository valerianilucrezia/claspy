library(ggplot2)
library(tidyverse)

results <- read.table('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/claspy/mean/results.tsv', sep = '\t', header = T)
data <- readRDS('~/Dropbox/projects/lr-cnc/data/smooth_data.RDS')

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
    plt <- plot_claspy(data, tmp, type = 'mean')  
  }
}



plot_claspy <- function(data, bp, true_bp, type = 'mean'){
  name <- paste0(tmp$mode, '_', tmp$thr, '_', tmp$ws)
  
  pred_res <- parse_res(tmp$res)
  dr <- paste0(type, '_dr')
  baf <- paste0(type, '_baf')
  
  
  plt <- ggplot() +
    geom_point(aes(x = data$order, y = data[[dr]]), size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylim(0,6) +
    ylab(dr) +
    xlab('pos') +
    
    ggplot() +
    geom_point(aes(x = data$order, y = data[[baf]]),size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylim(0,1) +
    ylab(baf) +
    xlab('pos') +
    
    ggplot() +
    geom_point(aes(x = data$order, y = data$gt_AF),size = 0.3, alpha = 0.2) +
    #geom_vline(aes(xintercept = true_bp), color = 'forestgreen', size = .8, alpha = 0.4, linetype = 'longdash') +
    geom_vline(aes(xintercept = chr_start_end), color = 'gray20', linetype = 'dotted') +
    geom_vline(aes(xintercept = pred_res), color = 'steelblue3') +
    scale_x_continuous(
      breaks = c(0, chr$half),
      labels = c("", gsub(pattern = 'chr', replacement = '', chr$CHROM))) + 
    ylab('VAF') +
    xlab('pos') +
    ylim(0,1) +
    plot_layout(nrow = 3, guides = 'collect') + plot_annotation(title = name) & theme_bw()
  
  ggsave(paste0('~/Desktop/orfeo_LTS/CNA/segmentation/res_nanopore/claspy/mean/claspy_',type, '_', name,'.png'), dpi = 400, plot = plt)
}



