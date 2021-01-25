
# somehow max.print default setting never works
options(max.print=200); # note this is not the maximum "rows" but maximum records
options(stringsAsFactors = FALSE)
## source('~/wzlib/Rutils/wzcore.R')
## source('~/wzlib/Rutils/wzplotlib.R')

# .First <- function(){
# }

load_default <- function() {
  library(grDevices)
  library(tidyverse)
  library(wheatmap)
  library(sesame)
  library(ggplot2)
  library(reshape2)
  library(devtools)
  library(readxl)
}

## .Last <- function(){
##   if(!any(commandArgs()=='--no-readline') && interactive()) {
##     hist_file <- Sys.getenv("R_HISTFILE")
##     if(hist_file=="") hist_file <- "~/Rhistory"
##     try({
##       timestamp(,prefix=paste("##--- SESSION_END --- [",getwd(),"] ",sep=""))
##       tmpfile <- tempfile("Rrawhist")
##       savehistory(tmpfile)
##       file.append(hist_file, tmpfile) 
##       unlink(tmpfile) 
##     })
##   }
## }
