get_true_bp_list <- function(real_bp){
  true_bp = list()
  for (s in real_bp$sim_id){
    tmp <- real_bp %>% filter(sim_id == s) %>% pull(real_bp)
    tmp <- strsplit(tmp, ', ') %>% unlist() 
    first <- strsplit(tmp[1], "\\[")[[1]][2]
    last <- strsplit(tmp[length(tmp)], "\\]")[[1]][1]
    tmp <- replace(tmp, 1, first)
    tmp <- replace(tmp, length(tmp), last)
    
    new_tmp <- c()
    for (idx in seq(1, length(tmp)-1)){
      if (as.integer(tmp[idx]) == (as.integer(tmp[idx+1]) -1)){
        c <- NULL
      } else {
        new_tmp <- c(new_tmp, as.integer(tmp[idx]))
      }
    }
    new_tmp <- c(new_tmp, as.integer(tmp[length(tmp)]))
    true_bp[[s]]$bp <- new_tmp
    
    
    len <- c()
    for (idx in seq(1, length(new_tmp)-1)){
      tmp <- as.integer(new_tmp[idx+1]) - as.integer(new_tmp[idx])
      len <- c(tmp, len)
    }
    true_bp[[s]]$len <- len
  }
  return(true_bp)
}


parse_res <- function(tmp){
  tmp <- strsplit(tmp, ', ') %>% unlist() 
  first <- strsplit(tmp[1], "\\[")[[1]][2]
  last <- strsplit(tmp[length(tmp)], "\\]")[[1]][1]
  tmp <- replace(tmp, 1, first)
  tmp <- replace(tmp, length(tmp), last)
  
  
  new_tmp <- c()
  for (idx in seq(1, length(tmp)-1)){
    if (as.integer(tmp[idx]) == (as.integer(tmp[idx+1]) -1)){
      c <- NULL
    } else {
      new_tmp <- c(new_tmp, as.integer(tmp[idx]))
    }
  }
  new_tmp <- c(new_tmp, as.integer(tmp[length(tmp)]))
  return(new_tmp)
}
