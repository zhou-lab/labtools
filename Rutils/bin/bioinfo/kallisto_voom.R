#!/usr/bin/env r
## this is memory heavy, since it loads the bam into memory, allow 10G
## runEdgeR.r -g mm10 -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o 
suppressMessages(library(docopt))

"Usage:
  kallisto_voom.r -g REFVERSION [-G GTF] -a CONDITION1 -b CONDITION2 -A SNAMES1 -B SNAMES2 -o OUTPUT

Options:
  -g REFVERSION   reference version (e.g, hg19, mm10)
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)   
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A SNAMES1      sample names for condition1 (comma-separated, kallisto/sname/abundance.tsv)
  -B SNAMES2      sample names for condition2 (comma-separated, kallisto/sname/abundance.tsv)
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
sink(stderr())
## cat(str(opt))

cat(format(Sys.time()), "Read samples...\n")
sample_list_n = c("PL_female_mut", "PL_male_mut", "PL_N9_female_mut", "PL_N9_male_mut")
for (i in 1:length(sample_list_n)) {
    tmp = read.table(file=paste0("kallisto/",sample_list_n[i],"/abundance.tsv"), header=T)
    assign(sample_list_n[i], tmp)       # attach each data frame to the current environment
}

sample_list = mget(sample_list_n)

## give the list unique names 
sample_list_uni = Map(function(x, i) setNames(x, ifelse(names(x) %in% "target_id",
      names(x), sprintf('%s.%d', names(x), i))), sample_list, seq_along(sample_list))

full_kalli = Reduce(function(...) merge(..., by = "target_id", all=T), sample_list_uni)
 
tpm_vals = full_kalli[, grep("tpm", names(full_kalli))]
rownames(tpm_vals) = full_kalli$target_id

## then make a contrast matrix
groups = c("normal", "colon_cancer", "normal", "colon_cancer")
condition &lt;- model.matrix(~0 + groups)
colnames(condition) = c("normal", "colon_cancer")
cont_matrix = makeContrasts(norm_canc = normal - colon_cancer, levels = condition)

## this spits out an EList object 
v = voom(counts = tpm_vals, design = condition)
fit = lmFit(v, condition)
fit = contrasts.fit(fit, cont_matrix)
fit = eBayes(fit)
top_table = topTable(fit, n = 10000, sort.by = "p")

browser()

## TODO


cat(format(Sys.time()), "Retrieving annotations...\n")
suppressPackageStartupMessages(library(GenomicFeatures))
if (is.null(opt$G)) {
    txdb <- suppressWarnings(makeTxDbFromUCSC(genome=opt$g,tablename='ensGene'))
    gtf <- NULL
} else {
    txdb <- makeTxDbFromGFF(opt$G, format="gtf",circ_seqs=character())
    seqlevelsStyle(txdb) <- "UCSC"
    suppressPackageStartupMessages(library(rtracklayer))
    gtf <- import(opt$G)
}

## txdb is TxDb class from GenomicFeatures package
## txByGene is GRangesList
txByGene <- transcriptsBy(txdb,'gene')

## NOTE: edgeR needs raw counts, not FPKM
cat(format(Sys.time()), "Making table of counts...\n")
suppressPackageStartupMessages(library(GenomicAlignments))
toc1 <- data.frame(sapply(unlist(strsplit(opt$A,',')), function(bam.path) {
    cat("Loading", bam.path, "..."); flush.console()
    gal <- readGAlignments(bam.path)
    reads <- GRanges(seqnames = rname(gal), ranges=IRanges(start=start(gal),end=end(gal)), strand=strand(gal))
    ## unstranded library
    ## reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),end=end(reads)), strand=rep("*",length(reads)))
    counts <- countOverlaps(txByGene, reads)
    cat("OK\n"); flush.console()
    counts
}, USE.NAMES=T), stringsAsFactors=F)

toc2 <- data.frame(sapply(unlist(strsplit(opt$B,',')), function(bam.path) {
    cat("Loading", bam.path, "..."); flush.console()
    gal <- readGAlignments(bam.path)    # GAlignments objects
    reads <- GRanges(seqnames = rname(gal), ranges=IRanges(start=start(gal),end=end(gal)), strand=strand(gal))
    ## unstranded library
    ## reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),end=end(reads)), strand=rep("*",length(reads)))
    counts <- countOverlaps(txByGene, reads)
    cat("OK\n"); flush.console()
    counts
}, USE.NAMES=T), stringsAsFactors=F)

toc <- cbind(toc1, toc2)

cat(format(Sys.time()), "Normalizing ...")
suppressPackageStartupMessages(library(edgeR))
normFactors = calcNormFactors(as.matrix(toc))
cat("OK\n")

cat(format(Sys.time()), "Statistical testing...\n")
group <- rep(c(opt$a, opt$b), c(dim(toc1)[2],dim(toc2)[2]))
DGE <- DGEList(toc, lib.size=normFactors*colSums(toc), group=group)
if (dim(toc1)[2]>1 && dim(toc2)[2]>1) { # with replica
    DGE.withdisp=estimateCommonDisp(DGE)
} else {                                # no replica
    cat("\nWarning: no replica, use GLM to estimate dispersion...")
    DGE.withdisp=estimateGLMCommonDisp(DGE,robust=T,subset=NULL,method="deviance")
}
tested=exactTest(DGE.withdisp)

cat(format(Sys.time()), "Outputing ...\n")
if (is.null(gtf)) {
    ## Old way of retrieving annotation from OrgDb, this lacks non-coding element names
    if (opt$g == 'mm10') {
        suppressPackageStartupMessages(library(org.Mm.eg.db))
        annodb <- org.Mm.eg.db
    } else if (opt$g == 'hg19') {
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        annodb <- org.Hs.eg.db
    }
    symbol <- mapIds(annodb, keys=row.names(tested$table), column="SYMBOL",keytype="ENSEMBL",multiVals="first")
} else {
    meta <- unique(as.data.frame(elementMetadata(gtf)[,c("gene_id","gene_name")]))
    id2name <- new.env(parent=emptyenv())
    invisible(apply(meta,1,function(x) id2name[[x[1]]] <<- x[[2]]))
    symbol <- unlist(mget(row.names(tested$table), envir=id2name));
}
write.table(cbind(DGE.withdisp$counts, tested$table, data.frame(symbol=symbol)), file=opt$o,sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.")
sink()
