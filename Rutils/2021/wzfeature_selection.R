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
    if (is.null(grouping) && is.null(group2probes)) {
        stop('Please specify either grouping or group2samples')
    }
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


#' This is adapted from branchInferProbes, but better.
#' It looks for one-vs-rest methylation signatures
SEInferTissueSpecificProbes = function(se, branch,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    ## hyper=FALSE; scan_delta = FALSE; scan_auc = FALSE
    betas = assay(se)
    if(is.null(rownames(betas))) rownames(betas) <- seq_len(nrow(betas))
    branch_grouping = as.data.frame(colData(se))[[branch]]

    in_na = rowSums(is.na(betas[,branch_grouping==0])) / sum(branch_grouping==0)
    out_na = rowSums(is.na(betas[,branch_grouping==1])) / sum(branch_grouping==1)
    
    ## hyper
    ## in-group min
    m0 = apply(betas[,branch_grouping==0],1,
        function(xx) mean(head(sort(xx),n=5), na.rm=TRUE))
    ## out-group max
    m1 = apply(betas[,branch_grouping==1],1,
        function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE))
    delta_beta = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, b1);})
    dfHype = data.frame(delta_beta = delta_beta, auc = auc[names(delta_beta)],
        in_na = in_na[names(delta_beta)], out_na = out_na[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hyper")

    ## hypo
    ## in-group max
    m0 = apply(betas[,branch_grouping==0],1,
        function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE))
    ## out-group min
    m1 = apply(betas[,branch_grouping==1],1,
        function(xx) mean(head(sort(xx),n=5), na.rm=TRUE))
    delta_beta = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, b1);})
    dfHypo = data.frame(delta_beta = -delta_beta, auc = auc[names(delta_beta)],
        in_na = in_na[names(delta_beta)], out_na = out_na[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hypo")

    ## branchScanDelta = function(df, auc_max = 0.01) {
    ##     message(sprintf("Found %d NA in Hyper", sum(is.na(df$delta_beta))))
    ##     for (i in seq(0.1,0.9,by=0.05)) {
    ##         probes = df$Probe_ID[
    ##         (!is.na(df$delta_beta)) & df$delta_beta >= i & df$auc <= auc_max]
    ##         message(sprintf("delta: %f Found %d probes", i, length(probes)))
    ##     }
    ## }
    message(sprintf("Processed branch %s (N_CpGs=%d)", branch, nrow(betas)))
    ## branchScanDelta(dfHypo)
    ## branchScanDelta(dfHype)
    rbind(dfHypo, dfHype)
}

filterSignature <- function(se, n_max = 50) {
    sigs = metadata(se)$tissue_signatures
    ## remove duplicate
    sigs = sigs %>% group_by(Probe_ID) %>% top_n(1, delta_beta)
    sigs = sigs %>% group_by(branch, type) %>%
        arrange(desc(delta_beta)) %>% dplyr::filter(row_number() <= n_max)
    stopifnot(length(sigs$Probe_ID) == length(unique(sigs$Probe_ID)))
    se = se[sigs$Probe_ID,]
    rownames(sigs) = sigs$Probe_ID
    rowData(se) = sigs
    se
}

orderBranch <- function(se) {
    sigs = as.data.frame(rowData(se))
    stopifnot(length(sigs$Probe_ID) == length(unique(sigs$Probe_ID)))
    CellTypes = unique(colData(se)$CellType)
    sigs1 = sigs %>% dplyr::filter(type=="Hypo")
    a = do.call(cbind, lapply(split(colData(se)$Sample_ID, colData(se)$CellType), function(x) {
        rowMeans(assay(se)[,x]) }))
    a = do.call(cbind, lapply(split(sigs1$Probe_ID, sigs1$branch), function(x) {
        colMeans(a[x,,drop=FALSE], na.rm=T) }))
    a = a[CellTypes,]
    taken = c()
    for(i in seq_len(nrow(a))) {
        js = which(a[i,] < 0.5 & !(seq_len(ncol(a)) %in% taken))
        taken = c(taken, js[order(-colSums(a[,js,drop=FALSE] < 0.5))])
    }
    orderedBranch = colnames(a[,taken])
    sigs$branch = factor(sigs$branch, level=orderedBranch)
    sigs = sigs %>% arrange(type, branch)
    se = se[sigs$Probe_ID,]
    rownames(sigs) = sigs$Probe_ID
    rowData(se) = sigs
    se
}

colClusterTissueSE <- function(se) {
    ## column cluster
    cd = as.data.frame(colData(se))
    grouping = cd$CellType
    betas = assay(se)
    betas2 = do.call(cbind, lapply(unique(grouping), function(g) {
        if (sum(grouping == g) > 3) {
            column.cluster(betas[,grouping == g])$mat
        } else {
            betas[,grouping == g,drop=FALSE]
        }
    }))
    se[, colnames(betas2)]
}

rowClusterTissueSE <- function(se) {
    ## row cluster
    rd = as.data.frame(rowData(se))
    grouping = paste0(rd$branch,"_", rd$type)
    betas = assay(se)
    betas2 = do.call(rbind, lapply(unique(grouping), function(g) {
        if (sum(grouping == g) > 2) {
            row.cluster(betas[grouping == g,])$mat
        } else {
            betas[grouping == g,,drop=FALSE]
        }
    }))
    se[rownames(betas2),]
}
