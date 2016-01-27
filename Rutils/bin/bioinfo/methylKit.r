#!/usr/bin/env r
## pipeline methylKit
## e.g., bioinfo/methylKit.r summary ~/.Renv/versions/3.2.3/lib64/R/library/methylKit/extdata/test1.myCpG.txt ~/.Renv/versions/3.2.3/lib64/R/library/methylKit/extdata/control1.myCpG.txt

## How to get annotation? see https://github.com/al2na/methylKit
## You can download annotation files from UCSC table browser for your genome of interest. Go to [http://genome.ucsc.edu/cgi-bin/hgGateway]. On the top menu click on "tools" then "table browser". Select your "genome" of interest and "assembly" of interest from the drop down menus. Make sure you select the correct genome and assembly. Selecting wrong genome and/or assembly will return unintelligible results in downstream analysis.

## From here on you can either download gene annotation or CpG island annotation.

##   1.  For gene annotation, select "Genes and Gene prediction tracks" from the "group" drop-down menu. Following that, select "Refseq Genes" from the "track" drop-down menu. Select "BED- browser extensible data" for the "output format". Click "get output" and on the following page click "get BED" without changing any options. save the output as a text file.
##   2.  For CpG island annotation, select "Regulation" from the "group" drop-down menu. Following that, select "CpG islands" from the "track" drop-down menu. Select "BED- browser extensible data" for the "output format". Click "get output" and on the following page click "get BED" without changing any options. save the output as a text file.

suppressMessages(library(docopt))

"Usage:
  methylKit.r summary [-t TREATMENT] [-b SAMPLEID] <INPUTFILE>... [-o OUTPUT]
  methylKit.r cluster [-t TREATMENT] [-b SAMPLEID] <INPUTFILE>... [-o OUTPUT]
  methylKit.r diff [-G GTF] -t <TREATMENT> -b <SAMPLEID> <INPUTFILE>... [-o OUTPUT]

Options:
  INPUTFILE       input files, must be of same length as treatment string (if supplied)
  -t TREATMENT    treatment string (e.g., 1,0) separated by ','
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)   
  -b SAMPLEID     sample id (','-separated, default to basename of file)
  -o OUTPUT       output dir [default: .]
"-> doc

opt <- docopt(doc)

suppressMessages(library(methylKit))

cat(format(Sys.time()), "Retrieving annotations...\n")
suppressPackageStartupMessages(library(GenomicFeatures))
if (is.null(opt$G)) {
    txdb <- suppressWarnings(makeTxDbFromUCSC(genome=opt$g,tablename='ensGene'))
} else {
    txdb <- makeTxDbFromGFF(opt$G, format="gtf",circ_seqs=character())
    seqlevelsStyle(txdb) <- "UCSC"
}

intronsByTx <- intronsByTranscript(txdb)

## cat(str(opt))

file.list <- opt$INPUTFILE

if (is.null(opt$SAMPLEID)) {
    opt$SAMPLEID = as.list(basename(opt$INPUTFILE))
} else {
    opt$SAMPLEID = as.list(unlist(strsplit(opt$SAMPLEID,",")))
}

n <- length(opt$INPUTFILE)

## if missing treatment
if (is.null(opt$TREATMENT)) {
    opt$TREATMENT = rep(0,times=n)
} else {
    opt$TREATMENT <- as.integer(unlist(strsplit(opt$TREATMENT,",")))
}

cat("Loading input..")
meth.objs <- read(as.list(opt$INPUTFILE), sample.id=opt$SAMPLEID, assembly="hg19", treatment=opt$TREATMENT)

if (opt$summary) {
    suppressMessages(library(graphics))
    sink(paste0(opt$o, "/meth.summary"))
    for (meth.obj in meth.objs) {

        ## calculate summary
        cat(meth.obj@sample.id,"\n","===\n", sep="")
        getMethylationStats(meth.obj, plot=F, both.strands=F)

        ## plot methstats
        bitmap(paste0(opt$o, "/", meth.obj@sample.id, ".methstats.png"), res=500)
        getMethylationStats(meth.obj, plot=T, both.strands=F)
        dev.off()

        ## plot coverage stats
        bitmap(paste0(opt$o, "/", meth.obj@sample.id, ".covstats.png"), res=500)
        getCoverageStats(meth.obj, plot=T, both.strands=F)
        dev.off()
    }
    sink()
}

cat("OK\n")

if (opt$diff) {
    meth <- unite(meth.objs, destrand=T)
    myDiff <- calculateDiffMeth(meth)

    inputfile <- c("methylKit/APCminSmadh3_vs_APCminNormal/WGBS_APCminSmadh3.methylKit","methylKit/APCminSmadh3_vs_APCminNormal/WGBS_APCmin_normal.methylKit")
    meth.raw.objs <- read(as.list(inputfile), sample.id=as.list(c("WGBS_APCminSmadh3","WGBS_APCmin_normal")), assembly="mm10", treatment=c(1,0)) # methylRaw class
    meth.base.obj <- unite(meth.raw.objs, destrand=T) # methylBase class
    diff <- calculateDiffMeth(meth)
    write.table(diff, "mydiff.tsv", sep="\t", quote=F, row.names=F)

    diffMethPerChr(diff, plot=F, qvalue.cutoff=0.01,meth.cutoff=25)
}

if (opt$cluster) {
    cat(str(meth.objs))
    ## getMethylationStats(myobj[[2]],plot=T,both.strands=F)
    ## meth = unite(myobj, destrand=T)
    ## cat(str(meth))
}
cat("Done.\n")
