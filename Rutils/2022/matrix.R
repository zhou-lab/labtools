cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
    cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
        f_row, f_col))
    cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
    namtx = is.na(mtx)
    good_row = rowSums(namtx) <= ncol(mtx) * f_row
    good_col = colSums(namtx) <= nrow(mtx) * f_col
    cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
    mtx[good_row, good_col]
}

imputeRowMean <- function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    mtx
}

#' Get most variable probes
bSubMostVariable <- function(betas, n=2000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}
