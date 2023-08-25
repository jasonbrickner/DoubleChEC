# ChEC-seq_Peak_Finder
## Introduction
ChIP-seq (chromatin immunoprecipitation followed by sequencing) is commonly used to identify genome-wide protein-DNA interactions. However, ChIP-seq often gives a low yield, which is not ideal for quantitative outcomes. An alternative method to ChIP-seq is ChEC-seq (Chromatin endogenous cleavage with high-throughput sequencing). In this method, the endogenous TF (transcription factor) of interest is fused with MNase (micrococcal nuclease) that non-specifically cleaves DNA near binding sites. Compared to the [original ChEC-seq method](https://www.nature.com/articles/ncomms9733), the [modified version](https://sites.northwestern.edu/bricknerlab/) requires far less amplification. Since [MACS3](https://github.com/macs3-project/MACS/tree/master#introduction) failed to identify peaks in data generated from the modified ChEC-seq method, a new peak finder has been developed specifically for it.

There are three functions in the *`peak_finder/`*. `callpeaks()` is used to identify peaks from BAM files. `goanalysis()` is used to make GO (Gene Ontology) term plots from peaks. `bedtomeme()` is a wrapper function to perform [MEME analysis](https://meme-suite.org/meme/tools/meme) in R **after [MEME Suite](https://meme-suite.org/meme/doc/download.html) is installed locally**.  
## Usage
1. Clone the ChEC-seq_Peak_Finder Repo
```
git clone https://github.com/ChengzheDuan/ChEC-seq_Peak_Finder.git
```
2. In the *`peak_finder/`* folder, find and open the *`run program.R`* file.
3. Set the working directory to where all the files are located
4. Load the `callpeaks()` function with:
```
source("peak finder.R")
```
5. Beofre running the peak finder with your own BAM files, run:
```
callpeaks(folder_treatment = "files_treatment_rap1",
          folder_control = "files_control_rap1",
          outdir = "test_folder_nup2_nohis_foldchange1.2_0.0001"
)
```
to test whether the peak finder works correctly.

7. Once the function finished running without any errors, load the `goanalysis()` function with:
```
source("GO analysis.R")
```
and run:
```
goanalysis(bedfile = "./test_folder/peaks.bed",
           outdir = "test_folder"
           )
```
to test whether a GO term plot is produced.

8. Lastly, load the `bedtomeme()` wrapper function with:
```
source("run meme.R")
```
and run:
```
bedtomeme(bed2fasta_filepath = "/opt/local/libexec/meme-5.5.1/bed2fasta",
          fasta_output_filename = "test_folder/fasta4meme.fna",
          genome_file = "sacCer3.fna",
          bed_filepath = "test_folder/peaks.bed",
          meme_filepath = "/opt/local/bin/meme",
          meme_output_folder = "test_folder/meme_output")
```
to test whether MEME results are produced. **Since `bedtomeme()` is a wrapper function, you have to have MEME Suite installed locally before using it. See the bedtomeme section below for more details**
## callpeaks
Use `callpeaks` function to identify peaks from BAM files. 
## Questions?
If you have any questions, please contact xxx
