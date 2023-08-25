bedtomeme <- function(bed2fasta_filepath, bed2fasta_additional = "",
                      fasta_output_filename, genome_file, bed_filepath,
                      meme_filepath, meme_additional = "", meme_output_folder){
  
  bed2pastastring <- paste(bed2fasta_filepath, 
                           bed2fasta_additional,
                           paste0("-o ", fasta_output_filename),
                           bed_filepath,
                           genome_file
                           )
  
  system(bed2pastastring, intern = T)
  
  memestring <- paste(meme_filepath, fasta_output_filename,
                      "-dna", "-revcomp", "-nmotifs 3", meme_additional,
                      paste0("-o ", meme_output_folder)
                       )
  
  system(memestring, intern = T)
}

