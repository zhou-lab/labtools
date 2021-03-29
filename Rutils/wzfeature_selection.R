getSignatureUW <- function(
    betas, grouping, u_max = 0.2, m_min = 0.7,
    max_na_in = 0, max_na_out = 0) {
    
    groups <- unique(grouping)
    is_na <- is.na(betas)
    sigs <- lapply(groups, function(g) {
        m1 <- rowMeans(betas[,grouping==g], na.rm=T) < u_max
        m2 <- rowMeans(betas[,grouping!=g], na.rm=T) > m_min
        ps1 <- rowSums(is_na[,grouping==g]) <= max_na_in
        ps2 <- rowSums(is_na[,grouping!=g]) <= max_na_out
        names(which(m1 & m2 & ps1 & ps2)) })
    names(sigs) <- groups
    sigs
}

getSignatureUTopW <- function(
    betas, grouping, n=100,
    max_na_in = 0, max_na_out = 0) {
    
    groups <- unique(grouping)
    is_na <- is.na(betas)
    sigs <- lapply(groups, function(g) {
        mean1 <- rowMeans(betas[,grouping == g], na.rm=T)
        mean0 <- rowMeans(betas[,grouping != g], na.rm=T)
        ps1 <- rowSums(is_na[,grouping == g]) <= max_na_in
        ps2 <- rowSums(is_na[,grouping != g] <= max_na_out)
        head(names(sort((mean1 - mean0)[ps1 & ps2])), n=n)
    })
    names(sigs) <- groups
    sigs
}

#' Cluster beta value matrix with Signature
clusterWithSignatureW <- function(betas, grouping, sigs) {
    pbs <- do.call(c, lapply(names(sigs), function(g) {
        if (length(sigs[[g]]) > 5)
            rownames(row.cluster(betas[intersect(
                rownames(betas), sigs[[g]]),])$mat)
        else
            NULL
    }))
    spl <- do.call(c, lapply(names(sigs), function(g) {
        colnames(column.cluster(betas[,grouping == g])$mat)
    }))
    betas[pbs, spl]
}

#' cluster betas by sample grouping
#' Either grouping or group2samples needs to be specified
clusterWithColumnGroupingW <- function(betas, grouping=NULL, group2samples=NULL, ordered_groups = NULL) {
    if (is.null(grouping) && is.null(group2samples)) stop('Please specify either grouping or group2samples')
    if (is.null(group2samples)) {
        if (is.null(ordered_groups)) {
            ordered_groups = unique(grouping)
        }
        do.call(cbind, lapply(ordered_groups, function(g) {
            if (sum(grouping == g) > 3) {
                column.cluster(betas[,grouping == g])$mat
            } else {
                betas[,grouping == g,drop=FALSE]
            }
        }))
    } else {
        do.call(cbind, lapply(group2samples, function(x) {
            if (is.null(x)) {
                NULL;
            } else if (length(x) > 2) {
                column.cluster(betas[,x])$mat
            } else {
                betas[,x,drop=FALSE]
            }
        }))
    }
}

clusterWithRowGroupingW <- function(betas, grouping=NULL, group2probes = NULL) {
    if (is.null(grouping) && is.null(group2probes)) stop('Please specify either grouping or group2samples')
    if (is.null(group2probes)) {
        do.call(rbind, lapply(unique(grouping), function(g) {
            if (sum(grouping == g) > 2) {
                row.cluster(betas[grouping == g,])$mat
            } else {
                betas[grouping == g,,drop=FALSE]
            }
        }))
    } else {
        do.call(rbind, lapply(group2probes, function(x) {
            if (!is.null(x) && length(x) > 2) {
                row.cluster(betas[x,])$mat
            } else if (length(x) == 0) {
                NULL;
            } else {
                betas[x,,drop=FALSE]
            }
        }))
    }
}

addDummyRows = function(betas, group2probes) {
    nc = ncol(betas)
    mtx = do.call(rbind, lapply(group2probes, function(betas) rbind(betas[probes,], rep(NA, nc))))
    mtx[1:(nrow(mtx)-1),]
}
