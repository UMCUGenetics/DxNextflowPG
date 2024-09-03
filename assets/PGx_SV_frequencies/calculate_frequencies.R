library(ggplot2)
library(dplyr)


Calc_frequencies <- function(gene) {
  data <- read.csv(paste0("~/tools/DxNextflowPG/assets/PGx_SV_frequencies/",gene,".csv"))

  data |>
    group_by_at(gene) |>
    summarise(Count = n()) |>
    mutate(Freq = Count/n()) |>
    write.csv(paste0("~/tools/DxNextflowPG/assets/PGx_SV_frequencies/",gene,"_freqs.csv"))

}


Calc_frequencies("CYP2D6")
Calc_frequencies("CYP2B6")
