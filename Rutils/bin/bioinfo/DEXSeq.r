#!/usr/bin/env r

suppressPackageStartupMessages(library(docopt))

"Usage:
  DEXSeq.r [--singleEnd] -G GTF -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -G GTF          DEXSeq-specific GTF
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A BAMS1        path to bams for condition1 (comma-separated)
  -B BAMS2        path to bams for condition2 (comma-separated)
  --singleEnd     treat reads as single-end
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
# cat(str(opt))

bpaths1 <- unlist(strsplit(opt$A,','))
bpaths2 <- unlist(strsplit(opt$B,','))
cntfiles <- c(bpaths1, bpaths2)

sampleTable <- data.frame(
    row.names = basename(cntfiles),
    condition = rep(c(opt$a, opt$b), c(length(bpaths1), length(bpaths2))),
    libType = rep("paired-end", length(opt$A)+length(opt$B)))

suppressPackageStartupMessages(library(DEXSeq))
cat(format(Sys.time()), "Reading count tables...\n")
dxd <- DEXSeqDataSetFromHTSeq(
    cntfiles,
    sampleData = sampleTable,
    design = ~ sample + exon + condition:exon,
    flattenedfile = opt$G)

cat(format(Sys.time()), "Normalizing...\n")
dxd <- estimateSizeFactors(dxd)
cat(format(Sys.time()), "Estimating dispersion...")
dxd <- estimateDispersions(dxd)

cat(format(Sys.time()), "Testing DEU...\n")
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")

cat(format(Sys.time()), "Outputing...\n")
dxr1 <- DEXSeqResults(dxd)
write.table(dxr1, opt$o, sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.\n")
