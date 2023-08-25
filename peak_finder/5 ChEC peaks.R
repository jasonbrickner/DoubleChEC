chec_peaks <- function(enriched_peaks, w, lower_distance, upper_distance){
  
  # Find ChEC peaks (i.e., two peaks within a certain distance)
  chec_peaks <- data.frame(chr = character(),
                           position = integer(),
                           start_bp = integer(),
                           end_bp = integer())
  
  for(count in 1:(nrow(enriched_peaks) - 1)){
    if((enriched_peaks$position[count + 1] - enriched_peaks$position[count]) >= upper_distance ||
       (enriched_peaks$position[count + 1] - enriched_peaks$position[count]) <= lower_distance){
      
    }else{
      if((enriched_peaks$position[count + 1] - enriched_peaks$position[count] ) > 0){
        chr = enriched_peaks$chr[count]
        position = (enriched_peaks$position[count] + enriched_peaks$position[count + 1]) / 2
        start_bp = as.integer(enriched_peaks$position[count] - w/2) 
        end_bp = as.integer(enriched_peaks$position[count + 1] + w/2)
        chec_peak <- cbind(chr, position, start_bp, end_bp)
        
        chec_peaks <- rbind(chec_peaks, chec_peak)
      }
    }
  }
  
  return(chec_peaks)
}
