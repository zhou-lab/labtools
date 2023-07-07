reorderRowsByBranch <- function(se) {
    cd = colnames(colData(se))
    celltype_names = grep("^CellType", cd, value=T)
    branches = unique(rowData(se)$branch)
    branches = branches[order(sapply(strsplit(branches,"\\.in\\."), function(x) {
        min(sapply(celltype_names, function(celltype_name) {
            ## make.names allow for space in the column titles
            match(x[1], make.names(colData(se)[[celltype_name]])) }), na.rm=T)
    }))]
    se = rbind(do.call(rbind, lapply(branches, function(branch) {
        se[rowData(se)$branch==branch & rowData(se)$type=="Hyper",]
    })), do.call(rbind, lapply(branches, function(branch) {
        se[rowData(se)$branch==branch & rowData(se)$type=="Hypo",]
    })))
}

selectBackground <- function(se, ranked_bgs) {
    contrasts = do.call(c, lapply(unique(colData(se)$CellType), function(ct) {
        ct = make.names(ct)
        contrast = NULL
        n = 0
        for (bg in paste0(".in.", ranked_bgs)) {
            dn = sum(rowData(se)$branch == paste0(ct, bg))
            if(dn > 0) {
                n = n + dn
                contrast = c(contrast, paste0(ct, bg))
                if (n > 50) { break; }}}
        contrast
    }))
    ## always include .in.All for groups
    se[grepl(".in.All$", rowData(se)$branch) | rowData(se)$branch %in% contrasts,]
}

defineHierarchicalContrasts <- function(meta) {
    grouplevels = grep("CellType",colnames(meta), value=T)
    cmpList = do.call(c, lapply(seq_along(grouplevels), function(l1) {
        do.call(c, lapply(seq_along(grouplevels), function(l0) {
            gs0 = meta[grouplevels[l0]]
            gs1 = meta[grouplevels[l1]]
            gs0uniq = unique(gs0)
            gs0uniq = gs0uniq[gs0uniq != "NA"]
            gs1uniq = unique(gs1)
            gs1uniq = gs1uniq[gs1uniq != "NA"]
            do.call(c, lapply(gs0uniq, function(g0) {
                lapply(gs1uniq, function(g1) {
                    st = ifelse(gs0 == g0 & gs1 == g1, 0,
                         ifelse(gs0 == g0 & gs1 != g1 & gs1 != "NA", 1, 2))
                    ## flip the states so 0011 and 1100 can be deduplicated
                    if(sum(st!=2)>0 && st[st!=2][1] == 1) {
                        st1 = ifelse(st == 2, 2, 1-st)
                        st1 = paste0(st1, collapse="")
                    } else {
                        st1 = paste0(st, collapse="")
                    }
                    list(nm_in = g1, nm_all = g0, st = st, st1 = st1)
                })
            }))
        }))
    }))

    cmpList = cmpList[sapply(cmpList, function(x) {
        sum(x$st == 0) > 0 &&
            sum(x$st == 1) > 0 &&
            sum(x$st == 0) / (sum(x$st == 0)+sum(x$st == 1)) > 0.01
    })]

    ## deduplicate
    states = sapply(cmpList, function(x) x$st1)
    cmpList = lapply(split(cmpList, states), function(x) x[[1]])
    cmpList = cmpList[sapply(cmpList, function(x) x$nm_all)!="Other" &
                      sapply(cmpList, function(x) x$nm_in)!="Other"]

    branches = do.call(cbind, lapply(cmpList, function(x) x$st))
    colnames(branches) = sapply(cmpList, function(x) paste0(x$nm_in,".in.",x$nm_all))
    meta1 = cbind(meta[,c("Sample_ID", grep("CellType",colnames(meta), value=T))], branches)
    meta1
}

auc_wmw2 = function(labels, scores) {
    labels <- as.logical(labels)
    n1 <- sum(labels)
    n2 <- sum(!labels)
    R1 <- sum(rank(scores)[labels])
    U1 <- R1 - n1 * (n1 + 1)/2
    U1/(n1 * n2)
}

#' This is adapted from branchInferProbes, but better.
#' It looks for one-vs-rest methylation signatures
SEInferTissueSpecificProbes = function(se, branch,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    ## hyper=FALSE; scan_delta = FALSE; scan_auc = FALSE
    betas = assay(se)
    if(is.null(rownames(betas))) rownames(betas) <- seq_len(nrow(betas))
    branch_grouping = colData(se)[[branch]]

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


#' This is adapted from branchInferProbes, but better.
#' It looks for one-vs-rest methylation signatures
SEMUInferTissueSpecificProbes = function(se, branch,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    ## hyper=FALSE; scan_delta = FALSE; scan_auc = FALSE
    betas = apply(assay(se), 2, function(x) MU2betas(x, mincov=min(20,quantile(MU2cov(x), probs=0.5))))
    rownames(betas) <- rownames(se)
    cov = MU2cov(assay(se))
    if(is.null(rownames(betas))) rownames(betas) <- seq_len(nrow(betas))
    branch_grouping = colData(se)[[branch]]

    in_na = rowSums(is.na(betas[,branch_grouping==0, drop=FALSE])) /
        sum(branch_grouping==0)
    out_na = rowSums(is.na(betas[,branch_grouping==1, drop=FALSE])) /
        sum(branch_grouping==1)
    in_minCov = rowMins(cov[,branch_grouping==0,drop=FALSE])
    names(in_minCov) = rownames(cov)
    out_meanCov = rowMeans(cov[,branch_grouping==1,drop=FALSE])
    names(out_meanCov) = rownames(cov)
    
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
        in_minCov = in_minCov[names(delta_beta)],
        out_meanCov = out_meanCov[names(delta_beta)],
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
        in_minCov = in_minCov[names(delta_beta)],
        out_meanCov = out_meanCov[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hypo")

    rbind(dfHypo, dfHype)
}


filterSignatureSE <- function(se, n_max = 50) {
    ## remove duplicate and select top cgs
    sigs <- as.data.frame(rowData(se)) %>% mutate(Probe_Index=seq_len(nrow(se)))
    sigs = sigs %>% group_by(Probe_ID) %>% top_n(1, delta_beta) %>% sample_n(1)
    sigs = sigs %>% group_by(branch, type) %>% arrange(desc(delta_beta)) %>% dplyr::filter(row_number() <= n_max)
    stopifnot(length(sigs$Probe_ID) == length(unique(sigs$Probe_ID)))
    se1 <- se[sigs$Probe_Index,]
    se1
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
    ## se = se1
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
