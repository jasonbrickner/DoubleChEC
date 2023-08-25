run_DESeq2 <- function(files_treatment, files_control, local_maxima_df,
                       numcores, chromosomes,
                       w, stepsize, foldchange, p_value, outdir, uas_info, SGD_YeastGenes){
  
  
  
  # Prepare treatment df for DESeq2
  df_list <- list()
  for(ft in 1:length(files_treatment)){
    results <- mclapply(chromosomes, function(chromosomes){
      
      df_sw <- data.frame(chr = character(),
                          position = integer(),
                          Count = numeric())
      
      chr_sub <- subset(files_treatment[[ft]], chr == chromosomes)
      chr_ave <- data.frame(chr = chromosomes,
                            position = round(SlidingWindow("mean", chr_sub$site, window = w, step = stepsize), 0),
                            Count = SlidingWindow("sum", chr_sub$Count, window = w, step = stepsize))
      df_sw <- rbind(df_sw, chr_ave)
    }, mc.cores = numcores) %>% bind_rows()
    
    df_list[[ft]] <- results
  }
  
  
  df_treatment <- data.frame(chr = df_list[[1]]$chr,
                             position = df_list[[1]]$position,
                             Count_1 = df_list[[1]]$Count,
                             Count_2 = df_list[[2]]$Count)
  
  if(length(df_list) == 2){
    
  }else{
    for(i in 3:length(df_list)){
      count_col_name = paste("Count", i, sep = "_")
      df_treatment[, count_col_name] = df_list[[i]]$Count
    }
  }
  
  
  # Prepare control df for DESeq2
  df_control_list <- list()
  for(i in 1:length(files_control)){
    results_control <- mclapply(chromosomes, function(chromosomes){
      
      df_sw_control <- data.frame(chr = character(),
                                  position = integer(),
                                  Count = numeric())
      
      chr_sub_control <- subset(files_control[[i]], chr == chromosomes)
      chr_ave_control <- data.frame(chr = chromosomes,
                                    position = round(SlidingWindow("mean", chr_sub_control$site, window = w, step = stepsize), 0),
                                    Count = SlidingWindow("sum", chr_sub_control$Count, window = w, step = stepsize))
      df_sw_control <- rbind(df_sw_control, chr_ave_control)
    }, mc.cores = numcores) %>% bind_rows()
    
    
    df_control_list[[i]] <- results_control
    
  }
  
  
  df_control <- data.frame(chr = df_control_list[[1]]$chr,
                           position = df_control_list[[1]]$position,
                           Count_1 = df_control_list[[1]]$Count,
                           Count_2 = df_control_list[[2]]$Count)  
  
  if(length(df_control_list) == 2){
    
  }else{
    for(i in 3:length(df_control_list)){
      count_col_name = paste("Count", i, sep = "_")
      df_control[, count_col_name] = df_control_list[[i]]$Count
    }
  }
  
  
  # Merge df_treatment and df_control
  merged_count_df <- dplyr::inner_join(df_control, df_treatment,
                           by = c("chr", "position"),
                           suffix = c("_control", "_treatment")
                           )
  
  selected_peaks_df_with_merged_count_df <- inner_join(local_maxima_df, merged_count_df, 
                                                       by = c("chr", "position")) %>% 
    select(-mean_CPMn)
  
  # Prepare matrix input for DESeq2
  df_4_DESeq2 <- selected_peaks_df_with_merged_count_df %>% 
    unite(data = .,
          col = "united_position",
          chr, position,
          sep = "-")
  
  row.names(df_4_DESeq2) <- df_4_DESeq2$united_position
  df_4_DESeq2$united_position <- NULL
  
  colnames(df_4_DESeq2) <- c(paste("untreated", 1:length(df_list), sep = ""), 
                             paste("treated", 1:length(df_control_list), sep = ""))
  
  
  # Prepare sample information table input for DESeq2
  v1 <- c(paste("untreated", 1:length(df_list), sep = ""), 
          paste("treated", 1:length(df_control_list), sep = ""))
  v2 <- c(rep("untreated", length(df_control_list)), rep("treated", length(df_list)))
  coldata <- data.frame(sample = as.factor(v1),
                        condition = as.factor(v2))
  
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = df_4_DESeq2,
                                colData = coldata,
                                design = ~ condition)
  
  dds <- DESeq(dds)
  DDS_results <- results(dds, lfcThreshold = foldchange,
                         altHypothesis="greater",
                         contrast=c("condition","treated","untreated"))
  
  DDS_results_df <- as.data.frame(DDS_results)
  
  
  DDS_results_df_annotated <- DDS_results_df %>% 
    mutate(united_position = row.names(DDS_results_df))
  
  row.names(DDS_results_df_annotated) <- NULL
  DDS_results_df_annotated <- DDS_results_df_annotated %>% 
    separate(data = .,
             col = united_position,
             into = c("chr", "position"),
             sep = "-") %>% 
    mutate("chromStart" = as.integer(as.numeric(position) - w/2),
           "chromEnd" = as.integer(as.numeric(position) + w/2)) %>% 
    select(-"position")
  
  UAS_GR <- makeGRangesFromDataFrame(uas_info[],
                                     seqnames.field=c("chr"))
  
  W_peaks <- DDS_results_df_annotated %>% mutate(strand = "+")
  C_peaks <- DDS_results_df_annotated %>% mutate(strand = "-")
  
  DDS_results_df_annotated <- rbind(W_peaks, C_peaks) %>% as.data.table()
  
  DDS_results_df_annotated_GR <- makeGRangesFromDataFrame(DDS_results_df_annotated[],
                                                          seqnames.field=c("chr") )
  UAS_overlaps <- findOverlapPairs(unname(UAS_GR), unname(DDS_results_df_annotated_GR))
  
  UAS_peaks <- as.data.table(UAS_overlaps) %>% select(c(1:3, 5, 7, 8)) %>% 
    dplyr::rename(start = first.start) %>%
    dplyr::rename(end = first.end) %>%
    dplyr::rename(peak_start = second.start) %>%
    dplyr::rename(peak_end = second.end) %>%
    dplyr::rename(chr = first.seqnames) %>%
    dplyr::rename(strand = first.strand)
  
  UAS_peaks_merged <- data.table::merge.data.table(UAS_peaks, DDS_results_df_annotated,
                                                   by.x = c("chr", "peak_start", "peak_end", "strand"),
                                                   by.y = c("chr", "chromStart", "chromEnd", "strand"))
  
  UAS_peak_genes <- right_join(uas_info, UAS_peaks_merged, 
                               by = c("chr", "start", "end", "strand"), 
                               multiple = "all",
                               relationship = "many-to-many")
  
  colnames(UAS_peak_genes)[5] <- "Systematic_Name"
  
  # Import yeast gene information

  list_peaks <- data.table::merge.data.table(UAS_peak_genes, SGD_YeastGenes, 
                                             by = "Systematic_Name") %>% 
    select(-c("start", "end", "Primary_SGID")) %>% 
    select(c(2, 4:11, 3, 1, 12))
  

  DDS_results_df_annotated_edited <- DDS_results_df_annotated %>% 
    select(7, 8, 9)
  colnames(DDS_results_df_annotated_edited) <- c("chr", "peak_start", "peak_end")
  list_peaks_edited <- list_peaks %>% select(1, 2, 3)

  local_max_no_annotation <- dplyr::setdiff(DDS_results_df_annotated_edited, list_peaks_edited)
  
  local_max_no_annotation <- data.table::merge.data.table(local_max_no_annotation, DDS_results_df_annotated,
                                                          by.x = c("chr", "peak_start", "peak_end"),
                                                          by.y = c("chr", "chromStart", "chromEnd")) %>% 
    select(-"strand") %>% 
    unique()
  
  all_local_max <- full_join(list_peaks, local_max_no_annotation, 
                             by = join_by(chr, peak_start, peak_end, baseMean,
                                          log2FoldChange, lfcSE, stat,
                                          pvalue, padj))
  
  
  dir.create(outdir)
  
  write.xlsx(all_local_max, paste(outdir, "Peak Annotation.xlsx", sep = "/"), sheetName = "Local Maxima")
  
  # Filter DESeq2 results
  DDS_results_df_subset <- subset(DDS_results_df, DDS_results_df$padj < p_value)
  DDS_results_df_subset$united_position <- row.names(DDS_results_df_subset)
  row.names(DDS_results_df_subset) <- NULL
  DDS_results_df_final <- DDS_results_df_subset %>% separate(data = .,
                                                             col = united_position,
                                                             into = c("chr", "position"),
                                                             sep = "-") %>% 
    select(c("chr", "position"))
  DDS_results_df_final$position <- as.integer(DDS_results_df_final$position)
  
  
  # Add start and end base pairs to selected_peaks
  enriched_peaks <- DDS_results_df_final %>% 
    mutate(start_bp = as.integer(position - w/2),
           end_bp = as.integer(position + w/2))
  
  W_peaks_enriched_peaks <- enriched_peaks %>% mutate(strand = "+")
  C_peaks_enriched_peaks <- enriched_peaks %>% mutate(strand = "-")
  
  enriched_peaks_annotated <- rbind(W_peaks_enriched_peaks, C_peaks_enriched_peaks) %>% 
    as.data.table() %>% 
    select(-"position") %>% 
    dplyr::rename(ChromStart = start_bp,
                  Chromend = end_bp)
  enriched_peaks_annotated_GR <- makeGRangesFromDataFrame(enriched_peaks_annotated[],
                                                          seqnames.field=c("chr") )
  
  UAS_overlaps_enriched_peaks <- findOverlapPairs(unname(UAS_GR), unname(enriched_peaks_annotated_GR)) %>% 
    as.data.table()
  
  UAS_peaks_enriched_peaks <- as.data.table(UAS_overlaps_enriched_peaks) %>% select(c(1:3, 5, 7, 8)) %>% 
    dplyr::rename(start = first.start) %>%
    dplyr::rename(end = first.end) %>%
    dplyr::rename(peak_start = second.start) %>%
    dplyr::rename(peak_end = second.end) %>%
    dplyr::rename(chr = first.seqnames) %>%
    dplyr::rename(strand = first.strand)
  
  UAS_peaks_merged_enriched_peaks <- data.table::merge.data.table(UAS_peaks_enriched_peaks, DDS_results_df_annotated,
                                                   by.x = c("chr", "peak_start", "peak_end", "strand"),
                                                   by.y = c("chr", "chromStart", "chromEnd", "strand"))
  
  UAS_peak_genes_enriched_peaks <- right_join(uas_info, UAS_peaks_merged_enriched_peaks, 
                                              by = c("chr", "start", "end", "strand"), 
                                              multiple = "all",
                                              relationship = "many-to-many")
  
  colnames(UAS_peak_genes_enriched_peaks)[5] <- "Systematic_Name"
  

  list_peaks_enriched_peaks <- data.table::merge.data.table(UAS_peak_genes_enriched_peaks, SGD_YeastGenes, 
                                             by = "Systematic_Name") %>% 
    select(-c("start", "end", "Primary_SGID")) %>% 
    select(c(2, 4:11, 3, 1, 12))
  
  ##############
  UAS_peak_genes_enriched_peaks_edited <- list_peaks_enriched_peaks %>% 
    select(1, 2, 3)
  colnames(UAS_peak_genes_enriched_peaks_edited) <- c("chr", "peak_start", "peak_end")
  list_peaks_enriched_peaks_edited <- list_peaks_enriched_peaks %>% select(1, 2, 3)
  
  enriched_peaks_edited <- enriched_peaks %>% select(-2)
  colnames(enriched_peaks_edited) <- c("chr", "peak_start", "peak_end")
  enriched_peaks_no_annotation <- dplyr::setdiff(enriched_peaks_edited, UAS_peak_genes_enriched_peaks_edited)
  
  enriched_peaks_no_annotation <- data.table::merge.data.table(enriched_peaks_no_annotation, DDS_results_df_annotated,
                                                          by.x = c("chr", "peak_start", "peak_end"),
                                                          by.y = c("chr", "chromStart", "chromEnd")) %>% 
    select(-"strand") %>% 
    unique()
  
  all_enriched_peaks <- full_join(list_peaks_enriched_peaks, enriched_peaks_no_annotation, 
                                  by = join_by(chr, peak_start, peak_end, baseMean,
                                               log2FoldChange, lfcSE, stat,
                                               pvalue, padj))
  
  peakannotation_excel <- loadWorkbook(paste(outdir, "Peak Annotation.xlsx", sep = "/"))
  addWorksheet(peakannotation_excel, "Peaks after sMNase Filter")
  writeData(peakannotation_excel, "Peaks after sMNase Filter", all_enriched_peaks)
  saveWorkbook(peakannotation_excel, 
               paste(outdir, "Peak Annotation.xlsx", sep = "/"),
               overwrite = TRUE)
  
  return(enriched_peaks)
  
}