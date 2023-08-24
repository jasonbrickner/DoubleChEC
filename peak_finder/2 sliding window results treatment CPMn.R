swresults_treatment_cpmn <- function(files_treatment, chromosomes, numcores,
                                     w, stepsize){
  
  # Generate sliding window results of treatment count data 
  sliding_window_treatment_count <- data.frame(chr = character(),
                                               position = integer(),
                                               CPMn = numeric())
  
  treatment_df_1_2 <- data.table::merge.data.table(files_treatment[[1]],
                                                   files_treatment[[2]],
                                                   by = c("chr", "site"),
                                                   suffixes = c("_1", "_2")) %>% 
    select(-c("filename_1", "filename_2"))

  if(length(files_treatment) == 2){
    
    merged_df <- treatment_df_1_2 %>% 
      mutate(ave_CPMn = rowMeans(treatment_df_1_2 %>% select(., starts_with("CPMn")))) %>% 
      select(c("chr", "site", "ave_CPMn"))
    
    results <- mclapply(chromosomes, function(chromosomes){
      
      df_sw <- data.frame(chr = character(),
                          position = integer(),
                          CPMn = numeric())
      
      chr_sub <- subset(merged_df, chr == chromosomes)
      
      chr_ave <- data.frame(chr = chromosomes,
                            position = round(SlidingWindow("mean", chr_sub$site, window = w, step = stepsize), 0),
                            CPMn = SlidingWindow("sum", chr_sub$ave_CPMn, window = w, step = stepsize))
      df_sw <- rbind(df_sw, chr_ave)
    }, mc.cores = numcores) %>% bind_rows()
    
    setDT(results)
    
    sw_averages_treatment_cpmn <- results[, .(mean_CPMn = mean(CPMn)), by = .(chr, position)]
  }else{
    for(i in 3:length(files_treatment)){
      cpmn_col_name = paste("CPMn", i, sep = "_")
      treatment_df_1_2[, cpmn_col_name] = files_treatment[[i]]$CPMn
    }
    merged_df <- treatment_df_1_2 %>% 
      mutate(ave_CPMn = rowMeans(treatment_df_1_2 %>% select(., starts_with("CPMn")))) %>% 
      select(c("chr", "site", "ave_CPMn"))
    results <- mclapply(chromosomes, function(chromosomes){
      
      df_sw <- data.frame(chr = character(),
                          position = integer(),
                          CPMn = numeric())
      
      chr_sub <- subset(merged_df, chr == chromosomes)
      chr_ave <- data.frame(chr = chromosomes,
                            position = round(SlidingWindow("mean", chr_sub$site, window = w, step = stepsize), 0),
                            CPMn = SlidingWindow("sum", chr_sub$ave_CPMn, window = w, step = stepsize))
      df_sw <- rbind(df_sw, chr_ave)
    }, mc.cores = numcores) %>% bind_rows()
    
    setDT(results)
    
    sw_averages_treatment_cpmn <- results[, .(mean_CPMn = mean(CPMn)), by = .(chr, position)]
    
  }
  
  return(sw_averages_treatment_cpmn)
}
