cuts <- function(bam) {
  bam_name <- strsplit(bam, "/")[[1]][length(strsplit(bam, "/")[[1]])]
  full_sc_ch <- readRDS("empty yeast genome coordinates.rds")
  GA <- readGAlignments(bam)
  cuts <- data.frame(seqnames(GA), strand(GA), start(GA), end(GA), stringsAsFactors = FALSE)
  colnames(cuts) <- c("chr", "strand", "start", "end")
  cuts$site <- (NA)
  cuts[cuts$strand == "-", ]$site <- cuts[cuts$strand == "-", ]$end
  cuts[cuts$strand == "+", ]$site <- cuts[cuts$strand == "+", ]$start
  MR <- length(cuts$site)/1e6
  cuts <- data.table(cuts[,-c(2:4)], key = c("chr", "site"))
  df4 <- cuts[, list(.N), by = key(cuts)]
  rep_df <- dplyr::left_join(full_sc_ch, df4, by = c("chr", "site"))

  rep_df[is.na(rep_df)] <- 0
  rep_df$CPMn <- rep_df$N/MR
  rep_df$ChEC_scaled <- rep_df$CPMn/mean(rep_df$CPMn)
  colnames(rep_df)[3] <- "Count"
  rep_df <- rep_df[,-5]
  rep_df$filename <- bam_name
  print(rep_df)
}

