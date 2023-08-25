importuas <- function(){
  yeast_UASs <- readRDS("yeast upstream regions.rds")
  
  uas_info <- yeast_UASs %>%
    mutate(strand_direction = case_when(strand == "+" ~ TRUE,
                                        strand == "-" ~ FALSE))
  
  uas_info$start <- ifelse(uas_info$strand_direction, uas_info$up_start - 200, uas_info$up_start)
  uas_info$end <- ifelse(uas_info$strand_direction, uas_info$up_end, uas_info$up_end + 200)
  uas_info <- uas_info %>% select(1, 7:8, 2:3) %>% .[order(.$chr, .$start),]

  
  return(uas_info)
}