
library(dmrseq)
library(tidyverse)

B = as.matrix(read_tsv("CG_BS_B"))
Cov = as.matrix(read_tsv("CG_BS_C"))

M = round(B*Cov)
cgfile = read_tsv("~/references/mm10/")