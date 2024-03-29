tsneSE <- function(se, perplexity=30, seed=1) {
    library(Rtsne)
    set.seed(seed)
    se = cleanMatrixForClusterSE(se)
    mx = imputeRowMean(assay(se))
    ## samples = colnames(mx)
    tsne = Rtsne(t(mx), dims=2, perplexity=perplexity)
    df = as.data.frame(tsne$Y)
    colnames(df) = c("tSNE1", "tSNE2")
    df$sample = colnames(mx)
    cbind(df, as_tibble(colData(se))) #[samples,]))
}

pcaSE <- function(se) {
    se = cleanMatrixForClusterSE(se, f_row=0.7, f_col=0.7)
    mx = imputeRowMean(assay(se))
    samples = colnames(mx)
    pca = prcomp(t(mx))
    cbind(pca$x, as_tibble(colData(se)[samples,]))
    ## cbind(as_tibble(colData(se)), pca$x[samples,])
    ## cbind(as_tibble(colData(se)), pca$x[colData(se)$IDAT,]) # or this
}

umapSE <- function(se, seed=1) {
    library(umap)
    set.seed(seed)
    se = cleanMatrixForClusterSE(se, f_row=0.7, f_col=0.7)
    mx = imputeRowMean(assay(se))
    samples = colnames(mx)
    umapResult <- data.frame(umap(t(mx))$layout)
    colnames(umapResult) <- c("UMAP1", "UMAP2")
    cbind(umapResult, as_tibble(colData(se)[samples,]))
}

