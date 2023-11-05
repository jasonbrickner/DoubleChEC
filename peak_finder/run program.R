rm(list = ls())

# Set working directory
# Please make sure all the scripts and treatment/control folders are in the working directory
setwd("set working directory")

# Peak calling
source("peak finder.R")

callpeaks(folder_treatment = "files_treatment_rap1",
          folder_control = "files_control_rap1",
          outdir = "test_folder"
)

# GO term plot
source("GO analysis.R")

goanalysis(bedfile = "./test_folder/peaks.bed",
           outdir = "test_folder"
           )

# Generate motif logos from bed file 
# If you are using a version of MEME other than 5.5.3, update meme version
meme_version = "meme-5.5.3"
source("run meme.R")

bedtomeme(bed2fasta_filepath = paste0("/opt/local/libexec/", meme_version, "/bed2fasta"),
          # bed2fasta_additional = "-both",
          fasta_output_filename = "test_folder/fasta4meme.fna",
          genome_file = "sacCer3.fna",
          bed_filepath = "test_folder/peaks.bed",
          meme_filepath = "/opt/local/bin/meme",
          # meme_additional = "-markov_order 2",
          meme_output_folder = "test_folder/meme_output")







