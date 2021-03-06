
# somehow max.print default setting never works
options(max.print=200); # note this is not the maximum "rows" but maximum records
options(stringsAsFactors = FALSE)
## source('~/wzlib/Rutils/wzcore.R')
## source('~/wzlib/Rutils/wzplotlib.R')

# .First <- function(){
# }

load_default <- function() {
    suppressMessages(source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R'))
    suppressMessages(source('~/repo/wzlib/Rutils/wzcore.R'))
    suppressMessages(source('~/repo/wzlib/Rutils/wzplotlib.R'))
    suppressMessages(source('~/repo/wzlib/Rutils/wzseq.R'))
    suppressMessages(source('~/repo/wzlib/Rutils/wzfeature_selection.R'))
    suppressMessages(library(grDevices))
    suppressMessages(library(tidyverse))
    suppressMessages(library(wheatmap))
    suppressMessages(library(scales))
    suppressMessages(library(sesame))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(SummarizedExperiment))
    suppressMessages(library(ggplot2))
    suppressMessages(library(reshape2))
    suppressMessages(library(devtools))
    suppressMessages(library(readxl))
    suppressMessages(library(viridis))
    suppressMessages(library(ggsci))
    suppressMessages(library(circlize))
    suppressMessages(library(ggrepel))
    suppressMessages(library(egg))
}

ld <- load_default

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
