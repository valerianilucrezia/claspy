get_true_bp_list <- function(real_bp){
  true_bp = list()
  for (s in unique(real_bp$sim_id)){
    tmp1 <- real_bp %>% filter(sim_id == s)

    for (c in tmp1$comb){
      tmp2 <- tmp1 %>% filter(comb == c) %>% pull(real_bp)
      tmp2 <- strsplit(tmp2, ', ') %>% unlist() 
      first <- strsplit(tmp2[1], "\\[")[[1]][2]
      last <- strsplit(tmp2[length(tmp2)], "\\]")[[1]][1]
      tmp2 <- replace(tmp2, 1, first)
      tmp2 <- replace(tmp2, length(tmp2), last)
    
      new_tmp <- c()
      for (idx in seq(1, length(tmp2)-1)){
        if (as.integer(tmp2[idx]) == (as.integer(tmp2[idx+1]) -1)){
          x <- NULL
        } else {
          new_tmp <- c(new_tmp, as.integer(tmp2[idx]))
        }
      }
      new_tmp <- c(new_tmp, as.integer(tmp2[length(tmp2)]))
      true_bp[[s]][[c]]$bp <- new_tmp  
      
      len <- c()
      for (idx in seq(1, length(new_tmp)-1)){
        tmp <- as.integer(new_tmp[idx+1]) - as.integer(new_tmp[idx])
        len <- c(tmp, len)
      }
      true_bp[[s]][[c]]$len <- len
    }
  }
  return(true_bp)
}


parse_res <- function(tmp){
  # tmp <- tmp$res
  
  first <- strsplit(tmp, "\\[")[[1]][2]
  last <- strsplit(first, "\\]")[[1]][1]
  tmp <- strsplit(last, ', ') %>% unlist() 
  
  if (length(tmp) == 1){
    return(tmp)
    
  } else if(length(tmp) == 0){
    return(NULL)
    
  } else {
    
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
}
