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

MU2betas_adaptive <- function(MU, mincov=1, maxcov=20, probs=0.1) {
    ## select coverage threshold so that 1-probs CpGs have value
    betas <- apply(MU, 2, function(x) {
        MU2betas(x, mincov=max(mincov, min(maxcov,
            quantile(MU2cov(MU), probs=probs))))
    })
    rownames(betas) <- rownames(MU)
    betas
}
