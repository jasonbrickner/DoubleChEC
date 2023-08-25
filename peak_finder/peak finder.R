
callpeaks <- function(folder_treatment, folder_control, 
                      numcores = NULL,
                      w = 3, stepsize = w / 2 + 0.5,
                      foldchange = 1.7, p_value = 0.0001, 
                      lower_distance = 15, upper_distance = 50, 
                      outdir = "output",
                      localmax_output = F,
                      peaks_after_solmnase_filter = F){
  
  # Set working directory to the current working directory
  setwd(getwd())

  
  # Get treatment file paths/names
  FL_treatment <- list.files(folder_treatment, pattern = "\\.bam$")
  
  for(i in 1:length(FL_treatment)){
    
    FL_treatment[i] = paste(folder_treatment, FL_treatment[i], sep = "/")
  }

  # Get control file paths/names
  FL_control <- list.files(folder_control, pattern = "\\.bam$")
  
  for(i in 1:length(FL_control)){
    
    FL_control[i] = paste(folder_control, FL_control[i], sep = "/")
  }
  
  
  if((length(FL_treatment) < 2) | (length(FL_control) < 2)){
    
    return("At least duplicate files are required!")
    
  }else{
    
    # Load required libraries
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(plyr)
    library(tidyverse)
    library(data.table)
    library(parallel)
    library(rlist)
    library(stringr)
    library(evobiR)
    library(fs)
    library(readxl)
    library(matrixTests)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Sc.sgd.db)
    library(lpSolve)
    library(irr)
    library(DESeq2)
    library(ggpmisc)
    library(openxlsx)
  
    # Set the default number of cores to use to run the program 
    if(Sys.info()["sysname"] == "Windows"){
      numcores = 1
    }else{
      if(is.null(numcores)){
        numcores = detectCores() - 1
      }
    }

  
    full_sc_ch <- readRDS("empty yeast genome coordinates.rds")
    
    # Map the sites of cleavage using the first basepair from each read, scale over the average count per bp.
    # 1a cuts.R and 1b mapcuts.R are used to count the number of times each position 
    # in the genome appears as the first position of a sequence read
    source("./1a cuts.R")
    source("./1b mapcuts.R")

    
    # Import and trim treatment files
    files_treatment <- mapcuts(FL_treatment, numcores = numcores)
    
    message("Treatment files are imported")
    
    
    # Import and trim control files
    files_control <- mapcuts(FL_control, numcores = numcores)
    
    message("Control files are imported")
    

    # Get a list of chromosomes 
    chromosomes <- unique(files_treatment[[1]]$chr)

    
    # Generate sliding window results from treatment CPMn data
    # 2 sliding window results treatment CPMn.R is used to smooth and average 
    # the treatment files using the sliding window and stepsize set by the user
    source("./2 sliding window results treatment CPMn.R")
    
    sw_averages_treatment_cpmn <- 
      swresults_treatment_cpmn(files_treatment = files_treatment, 
                               chromosomes = chromosomes,
                               numcores = numcores, w = w, stepsize = stepsize
    )
    
    message("Sliding window results of treatment CPMn data are generated")

    
    # Select all local maxima in the smoothed treatment data
    source("./3 select local maxima.R")
    
    local_maxima_df <- 
      localmax(sw_averages_treatment_cpmn = sw_averages_treatment_cpmn)
    
    # For each window that is identified as a peak, run DESeq2 to identify if 
    # the window is statistically significant compared to Sol MNase
    source("./4 run DESeq2.R")
    
    # 6 UAS regions.R is used to convert the gene upstream500 info to upstream700
    source("./6 UAS regions.R")
    uas_info <- importuas()
    SGD_YeastGenes <- read.csv("GO SGD_YeastGenesList.csv") %>% as.data.table()
    enriched_peaks <- run_DESeq2(files_treatment = files_treatment, 
                                 files_control = files_control,
                                 local_maxima_df = local_maxima_df,
                                 w = w, stepsize = stepsize,
                                 foldchange = foldchange, p_value = p_value,
                                 numcores = numcores, chromosomes = chromosomes, 
                                 outdir = outdir,
                                 uas_info = uas_info, SGD_YeastGenes = SGD_YeastGenes)
    
    message("DESeq2 is finished running")

  
    # Find doublet peaks (two peaks between the lower_distance and upper_distance limits)
    source("./5 ChEC peaks.R")
    chec_peaks <- chec_peaks(enriched_peaks = enriched_peaks, 
                             lower_distance = lower_distance, 
                             upper_distance = upper_distance,
                             w = w)
    
    message("Peaks are identified")
    
    # Export local maxima df in .bed file if localmax_output if set as TRUE
    if(localmax_output){
      local_maxima_df_export <- local_maxima_df %>% select(-"mean_CPMn")
      if(w<8){
        local_maxima_df_export <- local_maxima_df_export %>% 
          mutate(start_bp = as.integer(position -5),
                 end_bp = as.integer(position + 5)) %>% 
          select(-"position") %>% 
          mutate(name = c(1:nrow(local_maxima_df_export)),
                 strand = ".")
      }else{
        local_maxima_df_export <- local_maxima_df_export %>% 
          mutate(start_bp = as.integer(position - w/2),
                 end_bp = as.integer(position + w/2)) %>% 
          select(-"position") %>% 
          mutate(name = c(1:nrow(local_maxima_df_export)),
                 strand = ".")
      }
      write.table(local_maxima_df_export, paste(".", outdir, "local maxima.bed", sep = "/"),
                  row.names = F, col.names = F, sep = "\t", quote = F)
    }
    

    
    # Export peaks filtered by sMNase if peaks_after_solmnase_filter is set as TRUE
    if(peaks_after_solmnase_filter){
      if(w<8){
        enriched_peaks_export <- enriched_peaks %>% 
          mutate(start_bp = as.integer(position -5),
                 end_bp = as.integer(position + 5)) %>% 
          select(-"position") %>% 
          mutate(name = c(1:nrow(enriched_peaks)),
                 strand = ".")
      }else{
        enriched_peaks_export <- enriched_peaks %>% 
          mutate(start_bp = as.integer(position - w/2),
                 end_bp = as.integer(position + w/2)) %>% 
          select(-"position") %>% 
          mutate(name = c(1:nrow(enriched_peaks)),
                 strand = ".")
      }
      write.table(enriched_peaks_export, paste(".", outdir, "peaks after sMNase filter.bed", sep = "/"),
                  row.names = F, col.names = F, sep = "\t", quote = F)
    }
    

    
    # Output: .bed file
    bed <- data.frame(chrom = chec_peaks$chr,
                      chromStart = as.integer(chec_peaks$start_bp),
                      chromEnd = as.integer(chec_peaks$end_bp),
                      name = c(1:length(chec_peaks$chr)))
    
    ChEC_peaks <- bed %>% mutate(strand = ".")
    
    # Annotate ChEC-seq peaks
    UAS_GR <- makeGRangesFromDataFrame(uas_info[],
                                       seqnames.field=c("chr"))
    
    W_peaks <- bed %>% mutate(strand = "+")
    C_peaks <- bed %>% mutate(strand = "-")
    
    peaks <- rbind(W_peaks, C_peaks) %>% as.data.table()
    
    peaks_GR <- makeGRangesFromDataFrame(peaks[],
                                         seqnames.field=c("chrom") )
    
    UAS_overlaps <- findOverlapPairs(unname(UAS_GR), unname(peaks_GR))
    
    UAS_peaks <- as.data.table(UAS_overlaps) %>% select(c(1:3, 5, 7, 8)) %>% 
      dplyr::rename(start = first.start) %>%
      dplyr::rename(end = first.end) %>%
      dplyr::rename(peak_start = second.start) %>%
      dplyr::rename(peak_end = second.end) %>%
      dplyr::rename(chr = first.seqnames) %>%
      dplyr::rename(strand = first.strand)
    
    UAS_peak_genes <- right_join(uas_info, UAS_peaks, 
                                 by = c("chr", "start", "end", "strand"), 
                                 multiple = "all",
                                 relationship = "many-to-many")
    
    colnames(UAS_peak_genes)[5] <- "Systematic_Name"
    

    
    list_peaks <- data.table::merge.data.table(UAS_peak_genes, SGD_YeastGenes, 
                                               by = "Systematic_Name") %>% 
      select(-c("start", "end", "Primary_SGID")) %>% 
      select(c(2:5, 1, 6))
    
    bed <- bed %>% select(-"name")
    colnames(bed) <- c("chr", "peak_start", "peak_end")
    UAS_peak_genes <- UAS_peak_genes %>% select(1, 6, 7)

    
    peaks_no_annotation <- dplyr::setdiff(bed, UAS_peak_genes)
    
    all_peaks <- full_join(list_peaks, peaks_no_annotation, by = join_by(chr, peak_start, peak_end))

    
    peakannotation_excel <- loadWorkbook(paste(outdir, "Peak Annotation.xlsx", sep = "/"))
    addWorksheet(peakannotation_excel, "ChEC-seq Peaks")
    writeData(peakannotation_excel, "ChEC-seq Peaks", all_peaks)
    saveWorkbook(peakannotation_excel, 
                 paste(outdir, "Peak Annotation.xlsx", sep = "/"),
                 overwrite = TRUE)
    
    
    write.table(ChEC_peaks, paste(".", outdir, "peaks.bed", sep = "/"),
                row.names = F, col.names = F, sep = "\t", quote = F)
    

    # Output: description file with parameter settings and filtering results
    text <- paste("Parameters:", "\nw:", w, "\nstepsize:", stepsize, 
                  "\nfoldchange:", foldchange,
                  "\np_value:", p_value,
                  "\nlower distance:", lower_distance,
                  "\nupper distance:", upper_distance,
                  "\n\nPeak filtering results:",
                  "\nNumber of local maxima were found", nrow(local_maxima_df),
                  "\nNumber of peaks were enriched compared to Soluble MNase",nrow(enriched_peaks), 
                  "\nNumber of ChEC-seq peaks were found", nrow(ChEC_peaks))
    
    
    writeLines(text, paste(".", outdir, "Statistics.txt", sep = "/"))
    
    
    message("Done!")
    
  }
}









