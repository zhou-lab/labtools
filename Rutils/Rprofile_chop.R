
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
    suppressMessages(library(parallel))
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
    ## suppressMessages(library(viridis))
    suppressMessages(library(ggsci))
    suppressMessages(library(parallel))
    # suppressMessages(library(circlize))
    # suppressMessages(library(ggrepel))
    # suppressMessages(library(egg))
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


## rowMax <- function(x) {apply(x,1,max)}

## suppressWarnings(library(ggplot2))
## suppressWarnings(library(reshape2))
## suppressWarnings(suppressPackageStartupMessages(library(readxl)))
## suppressWarnings(suppressPackageStartupMessages(library(devtools)))
## suppressWarnings(suppressPackageStartupMessages(library(wheatmap)))
## suppressWarnings(suppressPackageStartupMessages(library(sesame)))
## suppressWarnings(suppressPackageStartupMessages(library(dplyr, quiet=TRUE)))
## suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
## suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges, quiet=TRUE)))
## suppressWarnings(suppressPackageStartupMessages(library(limma)))
## theme_wz <- theme_set(theme_classic(15))
## theme_wz <- theme_update(
##   axis.line.x = element_line(colour = "black"),
##   axis.line.y = element_line(colour = "black"),
##   axis.text = element_text(colour='black'))

## theme_wz_light <- theme_set(theme_linedraw(15))
## theme_wz_light <- theme_update(
##   ## text = element_text(face='bold'),
##   ## axis.line.x = element_line(colour = "black", size=1.2),
##   ## axis.line.y = element_line(colour = "black", size=1.2),
##   axis.text = element_text(colour='black', size=15),
##   panel.border = element_rect(linetype = "solid", colour = "black", size=1.2),
##   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## theme_wz <- theme_set(theme_linedraw(16))
## theme_wz <- theme_update(
##   text = element_text(face='bold'),
##   ## axis.line.x = element_line(colour = "black", size=1.2),
##   ## axis.line.y = element_line(colour = "black", size=1.2),
##   axis.text = element_text(colour='black', size=16),
##   panel.border = element_rect(linetype = "solid", colour = "black", size=1.2),
##   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
