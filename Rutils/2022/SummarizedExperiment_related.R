cleanMatrixForClusterSE <- function(se, f_row = 0.5, f_col = 0.5) {
    mtx = assay(se)
    cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
        f_row, f_col))
    cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
    namtx = is.na(mtx)
    good_row = rowSums(namtx) <= ncol(mtx) * f_row
    good_col = colSums(namtx) <= nrow(mtx) * f_col
    cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
    se[good_row, good_col]
}

bSubMostVariableSE <- function(se, n=2000) {
    mtx = assay(se)
    std <- apply(mtx, 1, sd, na.rm=TRUE)
    se[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

bSubVariableSE <- function(se, min_beta_range=0.2) {
    mtx = assay(se)
    beta_ranges = rowMaxs(mtx, na.rm=T) - rowMins(mtx, na.rm=T)
    se[beta_ranges > min_beta_range,]
}

bSubAutosomeSE <- function(se) {
    autoprobes <- names(sesameData_getAutosomeProbes(
        inferPlatformFromProbeIDs(rownames(se))))
    se[rownames(se) %in% autoprobes,]
}

bothClusterSE <- function(se, nrow_max = 3000, ncol_max = 3000) {
    if (nrow(se) > nrow_max || ncol(se) > ncol_max) {
        stop(sprintf("Too many rows (%d, max: %d) or columns (%d, max: %d). Abort.",
            nrow(se), nrow_max, ncol(se), ncol_max))
    }
    ord = both.cluster(assay(se))
    se[ord$row.clust$order, ord$column.clust$order]
}

columnClusterSE <- function(se, ncol_max = 3000) {
    if (ncol(se) > ncol_max) {
        stop(sprintf("Too many columns (%d, max: %d). Abort.",
            ncol(se), ncol_max))
    }
    ord = column.cluster(assay(se))
    se[, ord$column.clust$order]
}

columnClusterWithinGroupSE <- function(se, vars=NULL) {
    if (is.null(vars)) { vars = c("ExperimentGroup","Tissue") }
    df = colData(se) %>% as_tibble
    se1 = do.call(cbind, lapply(split(seq_len(nrow(df)), df[,vars]), function(ind) {
        if (length(ind)<2) {
            se[,ind]
        } else {
            columnClusterSE(se[,ind])
        }
    }))
    se1        
}

rowClusterWithinGroupSE <- function(se, vars=NULL) {
    if (is.null(vars)) { vars = c("branch","type") }
    df = rowData(se) %>% as_tibble
    rowgroups = split(seq_len(nrow(df)), df[,vars])
    rowgroups = rowgroups[sapply(rowgroups, length)>0]
    se1 = do.call(rbind, lapply(rowgroups, function(ind) {
        if (length(ind)<2) {
            se[ind,]
        } else {
            rowClusterSE(se[ind,], nrow_max=5000)
        }
    }))
    se1        
}

rowClusterSE <- function(se, nrow_max = 3000) {
    if (nrow(se) > nrow_max) {
        stop(sprintf("Too many columns (%d, max: %d). Abort.",
            nrow(se), nrow_max))
    }
    ord = row.cluster(assay(se))
    se[ord$row.clust$order,]
}

buildSE <- function(betas, probedf, meta, tissue_color, branch_color) {

    ## bt = clusterWithRowGroupingW(betas, group2probes = branch2probes)
    branch2probes = split(probedf$Probe_ID, probedf$branch)
    ## row order is defined by branch_color
    branch_color = na.omit(branch_color)
    bt = clusterWithRowGroupingW(betas, group2probes = branch2probes[names(branch_color)])
    ## column order is defined by tissue_color
    bt = clusterWithColumnGroupingW(bt, grouping = meta$tissue, ordered_groups = names(tissue_color))
    rd = probedf[match(rownames(bt), probedf$Probe_ID),]
    rownames(rd) = NULL
    cd = meta[match(colnames(bt), meta$Sample_ID),]
    se = SummarizedExperiment(assays=list(betas=bt), rowData=rd, colData=cd)
    metadata(se)$tissue_color = tissue_color
    metadata(se)$branch_color = branch_color
    se
}

getSEfromChunks <- function(rds_fdr, Probe_IDs) {
    se = do.call(rbind, mclapply(list.files(rds_fdr), function(x) {
        se1 = readRDS(sprintf("%s/%s", rds_fdr, x))
        rownames(se1) = paste0(seqnames(se1),"_",start(se1))
        se1 = se1[rownames(se1) %in% Probe_IDs,]
    }, mc.cores=24))
    se[Probe_IDs,]
}
