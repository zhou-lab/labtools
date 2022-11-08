MU2betas <- function(MU, mincov = 5) {
    U = bitwAnd(MU, 0xffff)
    M = bitwShiftR(MU, 16)
    F = M / (M+U)
    F[is.nan(F)] = NA
    F[U+M < mincov] = NA
    dim(F) = dim(MU)
    dimnames(F) = dimnames(MU)
    F
}
