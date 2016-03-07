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
dxr <- DEXSeqResults(dxd)

dxr.sorted <- dxr[order(dxr$padj),]
write.table(dxr.sorted[,names(dxr.sorted)!="transcripts"], opt$o, sep="\t", quote=F, col.names=NA)

plot.dir <- file.path(dirname(opt$o), "plots")
dir.create(plot.dir, showWarnings = FALSE)
gnames <- unique(dxr[dxr$padj<0.01 & !is.na(dxr$padj), ]$groupID)
for (i in seq_along(gnames)) {
    if (i>=100) break
    gname <- gnames[i];
    cat(format(Sys.time()), "Plotting", gname, "...\n");
    pdf(paste0(plot.dir,"/DEU_",gsub("\\+","_",gname),".pdf"))
    plotDEXSeq(dxr, gname, legend=T, cex.axis=1.2, cex=1.3, lwd=2)
    dev.off()
}

cat(format(Sys.time()), "Done.\n")
