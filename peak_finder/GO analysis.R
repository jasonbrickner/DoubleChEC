goanalysis <- function(bedfile, outdir, fontsize = 14, termlength = 10){
  
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
  
  source("6 UAS regions.R")
  source("GO Figure 2 GO term analysis.R")
  
  #find overlaps between peaks and UASs
  uas_info <- importuas()
  
  UAS_GR <- makeGRangesFromDataFrame(uas_info[],
                                     seqnames.field=c("chr") )
  
  peaks <- fread(bedfile)[,-c(4:5)] %>% 
    dplyr::rename(chr = V1,
                  start = V2,
                  end = V3)
  W_peaks <- peaks %>% mutate(strand = "+")
  C_peaks <- peaks %>% mutate(strand = "-")
  
  peaks <- rbind(W_peaks, C_peaks)
  
  peaks_GR <- makeGRangesFromDataFrame(peaks[],
                                       seqnames.field=c("chr"))
  
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
  
  GO_plotter(hit_list = unique(UAS_peak_genes$gene), outdir = outdir, fontsize = fontsize, termlength = termlength)
  
}




