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
