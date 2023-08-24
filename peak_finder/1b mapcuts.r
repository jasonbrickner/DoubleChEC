mapcuts <- function(files, numcores){
  file_list <- files
  result <- mclapply(file_list, cuts, mc.cores = numcores)
  
  return(result)
}


