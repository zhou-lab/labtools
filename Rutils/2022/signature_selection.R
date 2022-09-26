#' This is adapted from branchInferProbes, but better.
#' It looks for one-vs-rest methylation signatures
SEInferTissueSpecificProbes = function(se, branch,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    ## hyper=FALSE; scan_delta = FALSE; scan_auc = FALSE
    betas = assay(se)
    if(is.null(rownames(betas))) rownames(betas) <- seq_len(nrow(betas))
    branch_grouping = as.data.frame(colData(se))[[branch]]

    in_na = rowSums(is.na(betas[,branch_grouping==0, drop=FALSE])) /
        sum(branch_grouping==0)
    out_na = rowSums(is.na(betas[,branch_grouping==1, drop=FALSE])) /
        sum(branch_grouping==1)
    
    ## hyper: in-group min
    m0 = apply(betas[,branch_grouping==0, drop=FALSE],1,
        function(xx) mean(head(sort(xx),n=5), na.rm=TRUE))
    ## out-group max
    m1 = apply(betas[,branch_grouping==1, drop=FALSE],1,
        function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE))
    delta_beta = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, b1);})
    dfHype = data.frame(delta_beta = delta_beta, auc = auc[names(delta_beta)],
        in_na = in_na[names(delta_beta)], out_na = out_na[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hyper")

    ## hypo: in-group max
    m0 = apply(betas[,branch_grouping==0, drop=FALSE],1,
        function(xx) mean(tail(sort(xx),n=5), na.rm=TRUE))
    ## out-group min
    m1 = apply(betas[,branch_grouping==1, drop=FALSE],1,
        function(xx) mean(head(sort(xx),n=5), na.rm=TRUE))
    delta_beta = m0 - m1
    auc = apply(betas, 1, function(b1) {
        br = branch_grouping[branch_grouping %in% c(0,1)];
        b1 = b1[branch_grouping %in% c(0,1)];
        auc_wmw2(br, b1);})
    dfHypo = data.frame(delta_beta = -delta_beta, auc = auc[names(delta_beta)],
        in_na = in_na[names(delta_beta)], out_na = out_na[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hypo")

    rbind(dfHypo, dfHype)
}

filterSignatureGR <- function(se, n_max = 50) {
    ## remove duplicate and select top cgs
    gr <- rowRanges(se)
    names(gr) <- NULL
    sigs = GR2bed(gr)
    sigs$Probe_ID = paste0(sigs$seqnames, "_", sigs$start) # to be deleted
    sigs$row_ID = seq_len(nrow(sigs))
    ## top_n returns all if tie, sample_n randomly selects 1
    sigs = sigs %>% group_by(Probe_ID) %>% top_n(1, delta_beta) %>% sample_n(1)
    sigs = sigs %>% group_by(branch, type) %>% 
        arrange(desc(delta_beta)) %>% dplyr::filter(row_number() <= n_max)
    stopifnot(length(sigs$Probe_ID) == length(unique(sigs$Probe_ID)))
    se1 <- sort(se[sigs$row_ID,])
    mcols(se1)$Probe_ID <- paste0(seqnames(se1),"_",start(se1)) # to be deleted
    se1
}

orderBranchGR <- function(se) {
    ## order probe blocks according to the cell type/column order
    ## se0 = se; se = se1
    sigs = as.data.frame(rowData(se))
    stopifnot(length(sigs$Probe_ID) == length(unique(sigs$Probe_ID)))
    CellTypes = unique(colData(se)$CellType)
    sigs1 = sigs %>% dplyr::filter(type=="Hypo")
    a = do.call(cbind, lapply(split(seq_len(ncol(se)), colData(se)$CellType), function(x) {
        rowMeans(assay(se)[,x,drop=FALSE]) }))
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
    ## column cluster samples within cell type
    ## se = se3
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
