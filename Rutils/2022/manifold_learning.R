tsneSE <- function(se, perplexity=30, seed=1) {
    library(Rtsne)
    set.seed(seed)
    mx = assay(se)
    mx = imputeRowMean(cleanMatrixForClusterW(mx))
    tsne = Rtsne(t(mx), dims=2, perplexity=perplexity)
    df = as.data.frame(tsne$Y)
    colnames(df) = c("tSNE1", "tSNE2")
    df$sample = colnames(mx)
    cbind(df, as_tibble(colData(se)))
}

pcaSE <- function(se) {
    pca = prcomp(t(assay(se)))
    cbind(as_tibble(colData(se)), pca$x[colData(se)$IDAT,])
}
