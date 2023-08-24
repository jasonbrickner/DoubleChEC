# ChEC-seq_Peak_Finder
## Introduction
ChIP-seq (chromatin immunoprecipitation followed by sequencing) is commonly used to identify genome-wide protein-DNA interactions. However, ChIP-seq often gives a low yield, which is not ideal for quantitative outcomes. An alternative method to ChIP-seq is ChEC-seq (Chromatin endogenous cleavage with high-throughput sequencing). In this method, the endogenous TF (transcription factor) of interest is fused with MNase (micrococcal nuclease) that non-specifically cleaves DNA near binding sites. Compared to the [original ChEC-seq method](https://www.nature.com/articles/ncomms9733), the [modified ChEC-seq method](https://sites.northwestern.edu/bricknerlab/) requires far less amplification. Since [MACS3](https://github.com/macs3-project/MACS/tree/master#introduction) failed to identify peaks from data generated from the modified ChEC-seq method, a new peak finder has been developed to identify peaks from ChEC-seq data. 

## callpeaks
Use `callpeaks` function to identify peaks from BAM files. 
## Questions?
If you have any questions, please contact xxx
