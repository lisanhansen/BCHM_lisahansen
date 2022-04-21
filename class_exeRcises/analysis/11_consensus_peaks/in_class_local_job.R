library(tidyverse)
library(GenomicRanges)

source("/scratch/Shares/rinnclass/CLASS_2022/JR/CLASS_2022/util/class_functions.R")
#source("/scratch/Shares/rinnclass/CLASS_2022/JR/CLASS_2022/util/class_2022_functions.R")
consensus_peaks <- create_consensus_peaks(broadpeakfilepath = "/scratch/Shares/rinnclass/CLASS_2022/data/peaks/") 

# export consensus peaks to CORRECT DIRECTORY!
# This takes a lot longer (~ 1.5 hrs)
for(i in 1:length(consensus_peaks)) {
  rtracklayer::export(consensus_peaks[[i]], 
                      paste0("/scratch/Shares/rinnclass/CLASS_2022/LisaH/BCHM_lisahansen/class_exeRcises/analysis/11_consensus_peaks/consensus_peaks_in_class/", 
                             names(consensus_peaks)[i], 
                             "_consensus_peaks.bed"))
}