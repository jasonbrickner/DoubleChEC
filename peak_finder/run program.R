rm(list = ls())

# Set working directory
# Please make sure all the scripts and treatment/control folders in the working directory
setwd("/Users/chengzhe/Library/CloudStorage/GoogleDrive-chengzheduan@gmail.com/My Drive/Lab_Project/peak_finder_in_functions_upstream700_final_version")

# Peak calling
source("peak finder.R")

# Enter the treatment folder name where treatment files are stored and control
# folder name where control files are stored
callpeaks(folder_treatment = "files_treatment_rap1_new",
          folder_control = "files_control_rap1_new",
          outdir = "test_folder"
          # p_value = 0.0001, foldchange = 1.2
)

# Bar plot of GO terms 
source("GO analysis.R")

# Set the filepath/filename of the bed file generated from callpeaks() above
# Set the plot output directory 
goanalysis(bedfile = "./test_folder/peaks.bed",
           outdir = "test_folder"
           # fontsize = 10
           # termlength = 10
           )


# Generate motif logos from bed file
source("run meme.R")
# Download the MEME Suite here: https://meme-suite.org/meme/doc/download.html
# Install the MEME Suite by following one of the methods here: https://meme-suite.org/meme/doc/install.html?man_type=web
# Quick Install using MacPorts is recommended for Mac users

# bed2fasta_filepath: local bed2fasta program filepath (it's usually hidden. Mac users
# can use shift + command + . to find it)
# fasta_output_filename: output fasta filepath/filename 
# genome_file: Yeast genome file name(it's usually sacCer3.fna)
# Download the yeast genome file here: https://meme-suite.org/meme/doc/download.html
# bed_filepath: filepath/filename of the bed file for motif discovery
# bed2fasta_additional: pass additional bed2fasta arguments listed here: https://meme-suite.org/meme/doc/bed2fasta.html

# meme_filepath: local MEME program filepath (It's usually hidden. Mac users
# can use shift + command + . to find it)
# meme_output_folder: name of the folder created to write MEME results
# meme_additional: pass additional MEME arguments listed here: https://meme-suite.org/meme/doc/meme.html

bedtomeme(bed2fasta_filepath = "/opt/local/libexec/meme-5.5.1/bed2fasta",
          # bed2fasta_additional = "-both",
          fasta_output_filename = "test_folder/fasta4meme.fna",
          genome_file = "sacCer3.fna",
          bed_filepath = "test_folder/peaks.bed",
          meme_filepath = "/opt/local/bin/meme",
          # meme_additional = "-markov_order 2",
          meme_output_folder = "test_folder/meme_output")






