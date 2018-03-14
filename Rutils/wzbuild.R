#!/usr/bin/env Rscript
library(R6)
target <- commandArgs(trailingOnly=TRUE)
format <- "html"
if (length(target) > 1) {
  format <- target[2]
  target <- target[1]
}

if (dir.exists(target)) {
  library(devtools)
  build_vignettes(target)
  document()
  check()
} else {
  if (format == 'pdf') {
    library(rmarkdown)
    render(target, output_format='pdf_document')
  } else {
    library(rmarkdown)
    render(target)
  }
}
