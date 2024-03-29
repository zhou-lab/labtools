mergeOv <- function(a, b) {
  index <- findOverlaps(a, b)
  ans <- a[queryHits(index)]
  mcols(ans) <- c(mcols(ans), mcols(b[subjectHits(index)]))
  ans
}

excludeOv <- function(a, b) {
  index <- findOverlaps(a, b)
  a[-queryHits(index)]
}

lojOv <- function(a, b) {
  index <- findOverlaps(a, b)
  ans <- logical(length(a))
  ans[queryHits(index)] <- TRUE
  ans
}

table2GR <- function(tbl, rangeOnly = FALSE) {
    colnames(tbl)[1] <- 'chrm'
    colnames(tbl)[2] <- 'beg'
    colnames(tbl)[3] <- 'end'
    chrms <- sort(unique(tbl$chrm))
    if(is.factor(chrms)) { chrms <- as.character(chrms) }
    gr <- GRanges(
        seqnames=tbl$chrm, ranges=IRanges(tbl$beg, tbl$end),
        seqinfo=Seqinfo(chrms))
    if (!rangeOnly) {
        mcols(gr) <- tbl[,4:ncol(tbl)]
    }
    sort(gr)
}

bed2GR <- function(bedfn, rangeOnly=FALSE) {
    bed <- read.table(bedfn, header=F, stringsAsFactors=F)
    gr <- table2GR(bed, rangeOnly=rangeOnly)
    start(gr) <- start(gr) + 1 # bed format is 0-based at start
    gr
}

GR2bed <- function(gr) {
    start(gr) <- start(gr)-1
    gr$probeID <- names(gr)
    as.data.frame(gr)
}
