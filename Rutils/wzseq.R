##   df <- data.frame(seqnames=seqnames(gr),
##                    starts=start(gr)-1,
##                    ends=end(gr),
##                    names=c(rep(".", length(gr))),
##                    scores=c(rep(".", length(gr))),
##                    strands=strand(gr))

setMethod("$", "GRanges", function(x, name) {
  elementMetadata(x)[, name] })
setMethod("$<-", "GRanges", function(x, name, value) {
  elementMetadata(x)[[ name ]] <- value
  return(x)
})

GR2bed <- function(gr, fnm, rangeOnly=FALSE, col.names=F) {
  start(gr) <- start(gr)-1
  if (rangeOnly) {
    mcols(gr) <- NULL
  }
  gr$probeID <- names(gr)
  write.table(as.data.frame(gr), file=fnm, quote=F, sep="\t", row.names=F, col.names=col.names)
}

## intersection
setGeneric('%i%', function(x, y) standardGeneric('%i%'))
setMethod('%i%', c('ANY','ANY'), function(x, y) {intersect(x, y)})
## union
setGeneric('%u%', function(x, y) standardGeneric('%u%'))
setMethod('%u%', c('ANY','ANY'), function(x, y) {union(x, y)})
## setdiff
setGeneric('%d%', function(x, y) standardGeneric('%d%'))
setMethod('%d%', c('ANY','ANY'), function(x, y) {setdiff(x, y)})
## subsetByOverlaps
setGeneric('%s%', function(x, y) standardGeneric('%s%') )
setMethod('%s%', c('GRanges','GRanges'), function(x, y) {subsetByOverlaps(x, y)})
## ( hESC.H3K4ME1 %s% hESC.P300 ) %u% ( hESC.H3K4ME3 %s% hESC.H3K27ME3)
## can be very handy at times, and all of these operators are rather naturally associative IMHO.


reverse.complement <- function(s) {
  chartr("ACGT","TGCA",s)
}
