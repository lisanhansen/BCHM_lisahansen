library(tidyverse)
library(GenomicRanges)

# source("/scratch/Shares/rinnclass/CLASS_2022/JR/CLASS_2022/util/intersect_functions.R")
source("/scratch/Shares/rinnclass/CLASS_2022/LisaH/BCHM_lisahansen/util/class_functions.R")

# run consensus peaks
consensus_peaks <- create_consensus_peaks("/scratch/Shares/rinnclass/CLASS_2022/data/peaks")

# export consensus peaks 
for(i in 1:length(consensus_peaks)) {
  rtracklayer::export(consensus_peaks[[i]], 
                      paste0("/scratch/Shares/rinnclass/CLASS_2022/LisaH/BCHM_lisahansen/class_exercises/analysis/consensus_peaks/", 
                             names(consensus_peaks)[i], 
                             "_consensus_peaks.bed"))
}
