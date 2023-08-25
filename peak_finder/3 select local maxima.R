localmax <- function(sw_averages_treatment_cpmn){
  
  CPMn_threshold <- mean(sw_averages_treatment_cpmn$mean_CPMn) * 3
  
  select_local_maxima <- sw_averages_treatment_cpmn %>%
    mutate("peaks" = ggpmisc:::find_peaks(mean_CPMn, ignore_threshold = 0, span = 3)) %>%
    mutate(true_peaks = case_when((mean_CPMn > CPMn_threshold & peaks == TRUE) ~ TRUE,
                                  (mean_CPMn < CPMn_threshold & peaks == FALSE) ~ FALSE,
                                  TRUE ~ NA)) %>%
    select(-peaks)
  
  # df with selected local maxima only 
  local_maxima_df <- select_local_maxima %>% filter(true_peaks == TRUE) %>% select(-true_peaks)
  
  return(local_maxima_df)
}