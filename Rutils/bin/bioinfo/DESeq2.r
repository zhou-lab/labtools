#!/usr/bin/env r
## this is memory heavy, since it loads the bam into memory, allow 10G at least

## == transformation ==
## rld <- rlog(dds, blind=FLASE)

## == interpretating the result ==
## The first column, baseMean, is a just the average of the normalized count values, dividing by size factors, taken over all samples in the DESeqDataSet. lfcSE is "log fold change Standard Error"

## == other ways of importing data ==
## ddsSE <- DESeqDataSet(summarizedExperiment, design=~cell+dex) # DESeqDataSet class
## ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
## fpkm()
suppressMessages(library(docopt))

"Usage:
  DESeq2.r [--singleEnd] [--ignoreStrand] -g REFVERSION [-G GTF] [-s] -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -g REFVERSION   reference version (e.g, hg19, mm10)
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)   
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A BAMS1        path to bams for condition1 (comma-separated)
  -B BAMS2        path to bams for condition2 (comma-separated)
  --singleEnd     treat reads as single-end
  --ignoreStrand  ignore strand (should not set in stranded experiment)
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
sink(stderr())
## cat(str(opt))

cat(format(Sys.time()), "Retrieving annotations...\n")
suppressPackageStartupMessages(library(GenomicFeatures))
if (is.null(opt$G)) {
    txdb <- suppressWarnings(makeTxDbFromUCSC(genome=opt$g,tablename='ensGene'))
} else {
    txdb <- makeTxDbFromGFF(opt$G, format="gtf",circ_seqs=character())
    seqlevelsStyle(txdb) <- "UCSC"
    suppressPackageStartupMessages(library(rtracklayer))
    gtf <- import(opt$G)
}

## txdb is TxDb class from GenomicFeatures package
## exonByGene is GRangesList
suppressPackageStartupMessages(library(GenomicAlignments))
exonByGene <- exonsBy(txdb, by='gene') # GRangeList

cat(format(Sys.time()), "Summarizing experiment and making table of counts...\n")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("Rsamtools"))
bam.paths1 <- unlist(strsplit(opt$A,','))
bam.paths2 <- unlist(strsplit(opt$B,','))
bamfiles <- BamFileList(c(bam.paths1, bam.paths2), yieldSize=2000000)

## summarizeOverlaps can also take GRanges converted from "readGAlignments()" and "GRanges()"
## e.g., gal <- readGAlignments("~/run1_H55TVBGXX/bam/WGBS_APCmin_tumor.bam")
## reads <- GRanges(seqnames=rname(gal), ranges=IRanges(start=start(gal),
##   end=end(gal)), strand=strand(gal))
## below uses BamFileList from Rsamtools
library("BiocParallel")
register(SerialParam())                 # run in serial
## paired-end, capture fragments, strand-specific
fragments <- T
if (opt$singleEnd) fragments <- F
se <- summarizeOverlaps(features=exonByGene, reads=bamfiles, mode=Union, singleEnd=opt$singleEnd, ignore.strand=opt$ignoreStrand, fragments=fragments) # RangedSummarizedExperiment class

cat(format(Sys.time()), "Statistical testing...\n")
group <- rep(c(opt$a, opt$b), c(length(bam.paths1), length(bam.paths2)))
colData(se) <- cbind(colData(se), DataFrame(cond=as.factor(group)))
se$cond <- relevel(se$cond, opt$a)
dds <- DESeqDataSet(se, design=~cond)   # DESeqDataSet class, can be "design = ~ cell + dex"
dds <- DESeq(dds)
res <- results(dds)

cat(format(Sys.time()), "Outputing ...\n")
if (is.null(gtf)) {
    if (opt$g == 'mm10') {
        suppressPackageStartupMessages(library(org.Mm.eg.db))
        annodb <- org.Mm.eg.db
    } else if (opt$g == 'hg19') {
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        annodb <- org.Hs.eg.db
    }
    symbol <- mapIds(annodb, keys=row.names(res), column="SYMBOL",keytype="ENSEMBL",multiVals="first")
} else {
    meta <- unique(as.data.frame(elementMetadata(gtf)[,c("gene_id","gene_name")]))
    id2name <- new.env(parent=emptyenv())
    invisible(apply(meta,1,function(x) id2name[[x[1]]] <<- x[[2]]))
    symbol <- unlist(mget(row.names(res), envir=id2name));
}
write.table(cbind(res, as.data.frame(fpkm(dds)), data.frame(symbol=symbol)), opt$o, sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.")
sink()
