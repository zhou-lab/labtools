
# somehow max.print default setting never works
options(max.print=200); # note this is not the maximum "rows" but maximum records
options(stringsAsFactors = FALSE)
source('~/wzlib/Rutils/wzcore.R')
source('~/wzlib/Rutils/wzplotlib.R')


Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/pandoc')

.First <- function(){
  if(!any(commandArgs()=='--no-readline') && interactive()) {
    ## ... startup code (elided)
    ## print start directory and timestamp for searchability
    library(utils)
    timestamp(,prefix=paste("##-- SESSION_START ---- [",getwd(),"] ",sep=""))
  }
  library(grDevices)
  pdf.options(width=3.5, height=4, onefile=FALSE)
}

.Last <- function(){
  if(!any(commandArgs()=='--no-readline') && interactive()) {
    hist_file <- Sys.getenv("R_HISTFILE")
    if(hist_file=="") hist_file <- "~/Rhistory"
    try({
      timestamp(,prefix=paste("##--- SESSION_END --- [",getwd(),"] ",sep=""))
      tmpfile <- tempfile("Rrawhist")
      savehistory(tmpfile)
      file.append(hist_file, tmpfile) 
      unlink(tmpfile) 
    })
  }
}
