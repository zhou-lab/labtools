upsample_smote <- function(label, mx, lvls, n=100, k=3) {
    library(FNN)
    syn = lapply(lvls, function(lvl) {
        mx1 = mx[,label == lvl]
        n1 = n-ncol(mx1)
        if (ncol(mx1) <= k+2) {
            k1 = ncol(mx1) - 2
        } else {
            k1 = k
        }
        if (n1 > 0) {
            query_k = cbind(sample.int(ncol(mx1),n1,replace=T), sample.int(k1, n1, replace=T))
            synthetic = apply(query_k, 1, function(qk) {
                query = mx1[,qk[1]]
                neighbors = get.knnx(data = t(mx1), query = t(query), k1 = k1 + 1)
                ## +1 skip first neighbor
                neighbor_index <- neighbors$nn.index[1, qk[2]+1]
                neighbor <- mx1[,neighbor_index]
                ## random weight
                lambda <- runif(n = nrow(mx1))
                query*(1-lambda) + neighbor*lambda
            })
            list(mx = synthetic, label = rep(lvl, n1))
        } else {
            list(mx = NULL, label = NULL)
        }
    })
    list(mx = cbind(mx, do.call(cbind, lapply(syn, function(x) x$mx))),
        label = c(label, do.call(c, lapply(syn, function(x) x$label))))
}

oversample <- function(df, col, nmax) {
    do.call(bind_rows, lapply(unique(df[[col]]), function(level) {
        subset_df = df %>% dplyr::filter(.data[[col]] == level)
        rows_to_add = nmax - nrow(subset_df)
        if (rows_to_add > 0) {
            subset_df = bind_rows(subset_df,
                subset_df %>% slice_sample(n = rows_to_add, replace = TRUE))
        }
        subset_df
    }))
}

