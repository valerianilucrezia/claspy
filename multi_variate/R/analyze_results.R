library(ggplot2)
library(patchwork)
library(tidyverse)

source('~/Desktop/orfeo_LTS/CNA/segmentation/claspy/multi_variate/R/utils.R')

results <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/results.tsv', sep = '\t') %>% as_tibble()
real_bp <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/bps.tsv', sep = '\t')
true_bp <- get_true_bp_list(real_bp)


all_true <- c()
all_false <- c()

acceptance_wd = 2
for (row in seq(1, nrow(results))){
  #row = 2125
  tmp <- results %>% slice(row)
  comb <- tmp$comb
  
  pred_res <- parse_res(tmp$res)
  true_res <- true_bp[[tmp$sim_id]][[comb]]$bp
  
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
                              pred_wrong = all_false) %>% 
  separate(comb, 
           sep = '_', 
           into = c('NA1', 'cov', 'NA2', 'pur'), 
           remove = F) %>% 
  select(-NA1, -NA2)
                              

colors_mode = c(mult = 'cadetblue4', max = 'coral3', sum = 'darkseagreen4')


results %>% 
  pivot_longer(cols =c(pred_corr, pred_wrong)) %>% 
  ggplot() +
  geom_boxplot(aes(x = name, y = value, fill = mode)) +
  scale_fill_manual(values =colors_mode) +
  facet_grid(mode ~ thr ~ ws, scales = 'free_y') +
  theme_bw() + 


results %>% 
  filter(mode != 'max') %>% 
  filter(ws != 50) %>% 
  filter(thr != 1e-5) %>% 
  pivot_longer(cols =c(pred_corr, pred_wrong)) %>% 
  ggplot() +
  geom_boxplot(aes(x = name, y = value, fill = mode)) +
  scale_fill_manual(values =colors_mode) +
  facet_grid(mode ~ thr ~ ws ~ sim_id, scales = 'free_y') +
  theme_bw() + 


results %>% 
  mutate(cov = as.integer(cov),
         pur = as.numeric(pur)) %>% 
  filter(ws != 50) %>% 
  filter(thr != 1e-5) %>% 
  filter(mode != 'max') %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_corr, x = as.factor(cov), fill = mode)) +
  scale_fill_manual(values =colors_mode) +
  xlab('coverage') + 
  facet_grid(mode ~ thr ~ ws ~ as.factor(pur))  +


results %>% 
  mutate(cov = as.integer(cov),
         pur = as.numeric(pur)) %>% 
  filter(ws != 50) %>% 
  filter(thr != 1e-5) %>% 
  filter(mode != 'max') %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_wrong, x = as.factor(cov), fill = mode)) +
  scale_fill_manual(values =colors_mode) +
  xlab('coverage') +
  facet_grid(mode ~ thr ~ ws ~ as.factor(pur))+ 
  plot_layout(nrow = 2, ncol = 2, guides = 'collect')
ggsave('~/Desktop/acc_2.png', dpi = 300, width = 20, height = 15, units = 'in')
