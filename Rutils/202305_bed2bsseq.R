library(bsseq)
library(data.table)
library(tidyverse)

bed_dir = commandArgs(trailingOnly=TRUE)[1]
snames = sub("_cg.bed.gz","",grep("_cg.bed.gz", list.files(bed_dir), value=TRUE))

cgfiles = lapply(snames, function(sname) {
    fread(sprintf("%s/%s_cg.bed.gz", bed_dir, sname))
})

B = do.call(cbind, lapply(cgfiles, function(x) as.numeric(x$V4)))
Cov = do.call(cbind, lapply(cgfiles, function(x) as.integer(x$V5)))
Cov[is.na(Cov)] = 0
M = round(B*Cov)
M[is.na(M)] = 0
cgfile = with(cgfiles[[1]], tibble(chrm=V1, beg=V2, end=V3))

# idx = complete.cases(M)
# M = M[idx,]
# Cov = Cov[idx,]
# cgfile = cgfile[idx,]

bs <- BSseq(chr = cgfile$chrm, pos = cgfile$beg+1, M = M, Cov = Cov, sampleNames = snames)
saveRDS(bs, file="tmp/dmr/BSseq.rds")

