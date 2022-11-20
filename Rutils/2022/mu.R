MU2betas <- function(MU, mincov = 5) {
    library(bitops)
    U = bitAnd(MU, 0xffff)
    M = bitShiftR(MU, 16)
    F = M / (M+U)
    F[is.nan(F)] = NA
    F[U+M < mincov] = NA
    dim(F) = dim(MU)
    dimnames(F) = dimnames(MU)
    F
}

MU2cov <- function(MU, maxCov = 10000) {
    library(bitops)
    Cov = bitAnd(MU, 0xffff) + bitShiftR(MU, 16)
    Cov = pmin(Cov, maxCov)
    dim(Cov) = dim(MU)
    dimnames(Cov) = dimnames(MU)
    Cov
}
