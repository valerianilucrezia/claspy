library(patchwork)
library(dplyr)
library(ggplot2)

source('~/Desktop/orfeo_LTS/CNA/segmentation/claspy/multi_variate/R/utils.R')

results <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/results.tsv', sep = '\t') %>% as_tibble()
real_bp <- read.csv('~/Desktop/orfeo_LTS/CNA/segmentation/res_races/claspy/bps.tsv', sep = '\t')
true_bp <- get_true_bp_list(real_bp)

path_Data <- '~/Desktop/orfeo_LTS/CNA/segmentation/sim_data_races/data_races/'
#cna_color <- c('lightsteelblue', 'rosybrown2', 'wheat2', 'snow3','darkseagreen')
# t = 1e-15
# w = 5
# m = 'sum'
# sim = 'sim_10'

for (t in c(1e-15, 1e-10, 1e-5)){  
  for(w in c(5,10,50)){
    for (m in c('sum', 'mult', 'max')){
      name = paste0(m, '_', w, '_', t)
      print(name)
      
      plots <- list()
      for (s in seq(1,20)){
        print(s)
        
        sim <- paste0('sim_', s)
        path <- paste0(path_Data, sim)
        combinations <- list.dirs(path, full.names = F)
        combinations <- combinations[combinations != '']

        plottini <- lapply(combinations, function(c) {
          print(c)
          new_path <- paste0(path, '/', c)
          smooth_data <- readRDS(paste0(new_path, '/smooth_snv.RDS')) %>% arrange(pos)
          
          tb <- true_bp[[sim]][[c]]$bp
          pred <- results %>% 
            filter(sim_id == sim, 
                   comb == c) %>% 
            filter(thr == t) %>% 
            filter(mode == m) %>% 
            filter(ws == w)
          
          if (nrow(pred) > 0){
            pred_bp <- parse_res(pred$res)  
            if (length(pred_bp) == 1){
              pred_bp <- as.numeric(pred_bp)
            }
            
            plottino <- ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_dr), size = 0.5) + 
              geom_vline(aes(xintercept = tb), color = 'forestgreen', size = 1.5, alpha = 0.3) +
              geom_vline(aes(xintercept = pred_bp), color = 'steelblue') +  
              ggtitle(paste0(sim,'   ', c)) +
              xlab('pos') + 
              ylab('DR') +
              
              ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_baf), size = 0.5) + 
              geom_vline(aes(xintercept = tb), color = 'forestgreen', size = 1.5, alpha = 0.3) +
              geom_vline(aes(xintercept = pred_bp), color = 'steelblue') +  
              xlab('pos') +
              ylab('BAF') +
              
              ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_baf), size = 0.5) + 
              geom_vline(aes(xintercept = tb), color = 'forestgreen', size = 1.5, alpha = 0.3) +
              geom_vline(aes(xintercept = pred_bp), color = 'steelblue') +  
              xlab('pos') +
              ylab('VAF')  + plot_layout(nrow = 3) & theme_bw() 
          
            } else {
            plottino  <- ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_dr), size = 0.5) + 
              ggtitle(paste0(sim,'   ', c)) +
              xlab('pos') + 
              ylab('DR') +
              
              ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_baf), size = 0.5) + 
              xlab('pos') +
              ylab('BAF') +
              
              ggplot() +
              geom_point(aes(x = seq(1, nrow(smooth_data)), y = smooth_data$mean_baf), size = 0.5) + 
              xlab('pos') +
              ylab('VAF')  + plot_layout(nrow = 3) & theme_bw() 
          }
        })
        plots[[sim]] <- plottini
      }
      for (s in seq(1,20)){
        sim <- paste0('sim_', s)
        plt <- wrap_plots(plots[[sim]], ncol = 3)
        dir.create(paste0('~/Desktop/sims_plot/', sim, '/'), recursive = T, showWarnings = F)
        ggsave(filename = paste0('~/Desktop/sims_plot/', sim, '/', name, '.png'), plot = plt, height = 20, width = 10, dpi = 300, units = 'in')
        
      }
      
    }
  }
}
