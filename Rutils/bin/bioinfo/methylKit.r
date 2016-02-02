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
  methylKit.r diff -t <TREATMENT> -b SAMPLEID <INPUTFILE>... [-o OUTPUT] [-w CPGIBED] [-g REFSEQGENE]

Options:
  INPUTFILE       input files, must be of same length as treatment string (if supplied)
  -t TREATMENT    treatment string (e.g., 1,0) separated by ','
  -b SAMPLEID     sample id (','-separated, default to basename of file)
  -w CPGIBED      the UCSC table for CPGIBED
  -g REFSEQGENE   the UCSC table for REFSEQGENE
  -o OUTPUT       output dir [default: .]
"-> doc

opt <- docopt(doc)

suppressMessages(library(methylKit))

## cat(str(opt))

## sample id
if (is.null(opt$b)) {
    sampleid <- as.list(basename(opt$INPUTFILE))
} else {
    sampleid <- as.list(unlist(strsplit(opt$b,",")))
}

n <- length(opt$INPUTFILE)

## treatment
if (is.null(opt$t)) {
    treatment <- rep(0,times=n)
} else {
    treatment <- as.integer(unlist(strsplit(opt$t,",")))
}

cat("Loading input..")
meth.raw.objs <- read(as.list(unlist(strsplit(opt$INPUTFILE,","))), sample.id=sampleid, assembly="NA", treatment=treatment) # methylRawList class

if (opt$summary) {
    suppressMessages(library(graphics))
    sink(paste0(opt$o, "/meth.summary"))
    for (meth.obj in meth.raw.objs) {

        ## calculate summary
        cat(meth.obj@sample.id,"\n","===\n", sep="")
        getMethylationStats(meth.obj, plot=F, both.strands=F)

        ## plot methstats
        pdf(paste0(opt$o, "/", meth.obj@sample.id, ".methstats.pdf"))
        getMethylationStats(meth.obj, plot=T, both.strands=F)
        dev.off()

        ## plot coverage stats
        pdf(paste0(opt$o, "/", meth.obj@sample.id, ".covstats.pdf"))
        getCoverageStats(meth.obj, plot=T, both.strands=F)
        dev.off()
    }
    sink()
}

cat("OK\n")

if (opt$diff) {
    ## meth.raw.objs <- read(as.list(c("methylKit/APCminSmadh3_vs_APCminNormal/WGBS_APCminSmadh3.methylKit","methylKit/APCminSmadh3_vs_APCminNormal/WGBS_APCmin_normal.methylKit")), sample.id=as.list(c("WGBS_APCminSmadh3","WGBS_APCmin_normal")), assembly="mm10", treatment=opt$TREATMENT) # methylRawList class
    cat("Merging..\n")
    meth.base.obj <- unite(meth.raw.objs, destrand=T) # methylBase class
    cat("Calculating differential expression..\n")
    diff <- calculateDiffMeth(meth.base.obj)          # methylDiff class
    cat("Outputing single-base level..\n")
    write.table(cbind(meth.base.obj, diff), paste0(opt$o, "/diffmeth.tsv"), sep="\t", quote=F, row.names=F)

    cat("Get hyper/hypo methylated bases")
    ## get hypermethylated bases
    diff25p.hyper <- get.methylDiff(diff,difference=25,qvalue=0.01,type="hyper")
    ## get hypo methylated bases
    diff25p.hypo <- get.methylDiff(diff,difference=25,qvalue=0.01,type="hypo")

    sink(paste0(opt$o, "/diff.summary"))
    ## cat("Differential methylation per chromosome:")
    ## diffMethPerChr(diff, plot=F, qvalue.cutoff=0.01,meth.cutoff=25)

    ## annotate genic parts
    if (!is.null(opt$g)) {
        gene.obj <- read.transcript.features(opt$g)
        diffAnn <- annotate.WithGenicParts(diff, gene.obj) # annotationByGenicParts class
        ## getAssociationWithTSS(diffAnn)
        
        cat("Percentage/number of differentially methylated regions that overlap with intron/exon/promoters:")
        getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
        cat("\n")
        
        ## The percentage of differentially methylated bases overlapping with exon/intron/promoters
        pdf(paste0(opt$o, "/diffmeth_annotation.pdf"))
        plotTargetAnnotation(diffAnn,precedence=TRUE, main="differential methylation annotation")
        dev.off()

        cat("Percentage of intron/exon/promoters that overlap with differentially methylated bases.")
        getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

        ## conglomerate promoter methylation
        ## promoters <- regionCounts(myobj,gene.obj$promoters)
        ## write.table(promoters, paste0(opt$o, "/promoter.meth.tsv"), quote=F, sep="\t", row.names=F)
    }
    sink()

    ## annotate cpg island
    if (!is.null(opt$w)) {
        cpg.obj <- read.feature.flank(opt$w, feature.flank.name=c("CpGi","shores")) # GenomicRangesList
        diffCpGann <- annotate.WithFeature.Flank(diff, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi", flank.name="shores") # annotationByFeature class
        print(diffCpGann)
        ## percentage of differentially methylated bases are on CpG islands, CpG island shores and other regions.
        pdf(paste0(opt$o, "/diff_CpG_methylation.pdf"))
        plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="percentage of diffmeth bases in each cat")
        dev.off()
    }
    cat("\n")
}

if (opt$cluster) {
    cat(str(meth.raw.objs))
    ## getMethylationStats(myobj[[2]],plot=T,both.strands=F)
    ## meth = unite(myobj, destrand=T)
    ## cat(str(meth))
}
cat("Done.\n")
