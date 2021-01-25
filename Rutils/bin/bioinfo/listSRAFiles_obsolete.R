#!/usr/bin/env Rscript
## SRAlistFiles.R SRP026197
library(SRAdb)
srafile = getSRAdbFile()

con = dbConnect('SQLite',srafile)
listSRAfile(commandArgs(trailingOnly=TRUE)[1],con)

