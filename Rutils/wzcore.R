## probeChromosome <- function(
##     probes,
##     platform=c('EPIC','HM450','MM285'),
##     refversion=c('hg19','hg38','mm10')) {

##     platform <- match.arg(platform)
##     refversion <- match.arg(refversion)
    
##     mft <- sesameDataGet(sprintf(
##         '%s.%s.manifest', platform, refversion))[probes]
##     as.character(GenomicRanges::seqnames(mft))
## }

## probeIsAuto <- function(
##     probes,
##     platform=c('EPIC','HM450','MM285'),
##     refversion=c('hg19','hg38','mm10')) {

##     platform <- match.arg(platform)
##     refversion <- match.arg(refversion)

##     mft <- sesameDataGet(sprintf(
##         '%s.%s.manifest', platform, refversion))[probes]
##     !(as.character(GenomicRanges::seqnames(mft)) %in% c('chrX','chrY'))
## }

printf <- function(...) cat(sprintf(...));

wz <- function(df) {
    df.dim <- dim(df);
    print(df[1:min(10,df.dim[1]),1:min(10,df.dim[2])]);
    printf("%d rows %d columns\n", df.dim[1], df.dim[2]);
}

pbsgen <- function(
  dest=NULL,
  commands=NULL, args=NULL,
  jobname='WandingJob',
  queue='default',
  ppn=5,
  memG=2,
  walltime=12,
  submit=FALSE) {
  
  sink(paste0(dest,'.pbs'))
  cat(sprintf("
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N %s
#PBS -e %s.pbs.stderr
#PBS -o %s.pbs.stdout
#PBS -q %s
#PBS -l nodes=1:ppn=%s,mem=%sgb,walltime=%s:00:00
Rscript -e \"load('%s.rda'); myfun(args);\"
", jobname, dest, dest, queue, ppn, memG, walltime, dest))
  sink()

  myfun <- commands
  args <- args
  save(myfun, args, file=paste0(dest,'.rda'))

  cat('pbsfile: ', dest, '.pbs\n', sep='')
  if(submit) {
    system(paste0('qsub ', dest, '.pbs'))
  }
}

wzbind.list <- function(x, cat.names=NULL) {
  if (!is.null(cat.names))
    names(x) <- cat.names
  if (is.null(names(x)))
    stop('list has no names, please provide cat.names=')
  data.frame(x=do.call(c, x), cat=rep(names(x), sapply(x,length)))
}

wzbind <- function(..., names=NULL) {
  x <- list(...)
  if (is.null(names))
    names <- paste0('V',1:length(x))
  data.frame(x=do.call(c, x), cat=rep(names, sapply(x, length)))
}

wzvenn <- function(..., dnames=NULL) {
  data <- list(...)
  if (is.null(dnames))
    dnames <- 1:length(data)
  names(data) <- dnames
  library(VennDiagram)
  g <- venn.diagram(data, filename = NULL)
  grid.newpage()
  pushViewport(viewport(unit(0.1,'npc'),unit(0.1,'npc'),unit(0.8,'npc'),unit(0.8,'npc'), just=c('left','bottom')))
  grid.draw(g)
}

mergeOv <- function(a, b) {
  index <- findOverlaps(a, b)
  ans <- a[queryHits(index)]
  mcols(ans) <- c(mcols(ans), mcols(b[subjectHits(index)]))
  ans
}

excludeOv <- function(a, b) {
  index <- findOverlaps(a, b)
  a[-queryHits(index)]
}

lojOv <- function(a, b) {
  index <- findOverlaps(a, b)
  ans <- logical(length(a))
  ans[queryHits(index)] <- TRUE
  ans
}

## http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
slideMean <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

smoothScatter1 <- function(...) {
    smoothScatter(..., nrpoints=0, colramp=colorRampPalette(c("white","white","lightblue","blue","green","yellow","orange","red","darkred"),space = "Lab"), col='blue')
}

tbk_pack_fromIDAT <- function(pfxs, out_dir, idx_dir='~/references/InfiniumArray', mc.cores=4) {
    source('https://raw.githubusercontent.com/zhou-lab/tbmate/master/scripts/tbmate.R')
    idx_dir=path.expand(idx_dir)
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    tmp <- mclapply(seq_along(pfxs), function(i) {
        sset <- readIDATpair(pfxs[i])
        cat(pfxs[i],sset@platform,'\n')
        betas <- getBetas(dyeBiasCorrTypeINorm(noob(sset)), mask=FALSE)
        pvals <- pval(sset)[names(betas)]
        NULL
    } , mc.cores=mc.cores)
}


#' Extract the probe type field from probe ID
#' This only works with the new probe ID system.
#' See https://github.com/zhou-lab/InfiniumAnnotation for illustration
#'
#' @param Probe_ID Probe ID
#' @return a vector of '1' and '2' suggesting Infinium-I and Infinium-II
#' probeID_designType("cg36609548_TC21")
probeID_designTypeW <- function(Probe_ID) {
    stopifnot(all(grepl('_', Probe_ID))) # make sure it's the new ID system
    vapply(Probe_ID, function(x) substr(
        strsplit(x,'_')[[1]][2],3,3), character(1))
}

isUniqProbeIDW <- function(Probe_ID) {
    all(grepl('_',Probe_ID))
}

extractDesignW <- function(design_str) {
    vapply(
        stringr::str_split(design_str, ','),
        function(x) stringr::str_split(x[[1]],';')[[1]][1], character(1))
}

getProbesByRegionW <- function(
    chrm, beg = 1, end = -1,
    platform = c('EPIC','HM450'),
    refversion = c('hg38','hg19')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    if (end < 0) {
        end <- sesameDataGet(paste0(
            'genomeInfo.', refversion))$seqInfo[chrm]@seqlengths
    }

    probes <- sesameDataGet(paste0(
        platform, '.probeInfo'))[[paste0('mapped.probes.', refversion)]]
    
    if (!(chrm %in% GenomicRanges::seqinfo(probes)@seqnames)) {
        stop('No probes found in this reference');
    }
    message(sprintf('Extracting probes from %s:%d-%d.\n', chrm, beg, end))
    target.region <- GenomicRanges::GRanges(chrm, IRanges::IRanges(beg, end))
    subsetByOverlaps(probes, target.region)
}

getProbesByChromosomeW <- function(
    chrms, platform = c('EPIC','HM450'),
    refversion=c('hg19','hg38')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    mft <- sesameDataGet(sprintf('%s.%s.manifest', platform, refversion))
    names(mft)[as.character(GenomicRanges::seqnames(mft)) %in% chrms]
}

getAutosomeProbesW <- function(
    platform=c('EPIC','HM450','MM285'),
    refversion=c('hg19','hg38','mm10')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    mft <- sesameDataGet(sprintf(
        '%s.%s.manifest', platform, refversion))
    names(mft)[!(as.character(
        GenomicRanges::seqnames(mft)) %in% c('chrX', 'chrY'))]
}

bSubMostVariableW <- function(betas, n=2000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

bSubProbesW <- function(betas, probes, exclude=FALSE) {
    if (is.null(dim(betas))) { # should also work for vector
        if (exclude) {
            betas[!(names(betas) %in% probes)]
        } else {
            betas[intersect(names(betas), probes)]
        }
    } else {
        if (exclude) {
            betas[!(rownames(betas) %in% probes),]
        } else {
            betas[intersect(rownames(betas), probes),]
        }
    }
}

bSubCompleteW <- function(betas) {
    if (is.null(dim(betas))) { # should also work for vector
        betas[!is.na(betas)]
    } else {
        betas[complete.cases(betas),]
    }
}

bSubNonNAW <- function(betas, frac=0.5) {
    if (is.null(dim(betas))) { # should also work for vector
        betas[sum(!is.na(betas)) > frac*length(betas)]
    } else {
        betas[rowSums(!is.na(betas)) > frac*ncol(betas),]
    }
}

bSubNonNAColW <- function(betas, frac=0.5) {
    if (is.null(dim(betas))) { # should also work for vector
        betas[sum(!is.na(betas)) > frac*length(betas)]
    } else {
        betas[,colSums(!is.na(betas)) > frac*nrow(betas)]
    }
}

bSubCpGW <- function(betas) {
    betas[grep('^cg', rownames(betas)),]
}

bSubNoMaskW <- function(betas, platform='EPIC', refversion='hg38') {
    platform = 'EPIC'
    refversion = 'hg38'
    mft <- sesameDataGet(sprintf('%s.%s.manifest', platform, refversion))
    p <- intersect(rownames(betas), names(mft[!mft$MASK_general]))
    betas[p,]
}

extraHasW <- function(sset, k) {
    k %in% names(extra(sset))
}

extraGetW <- function(sset, k) {
    extra(sset)[[k]]
}

extraSetW <- function(sset, k, v) {
    extra(sset)[[k]] <- v
    sset
}

cbind_betas_onCommonW <- function(...) {
    input <- list(...)
    common <- Reduce(intersect, lapply(input, rownames))
    do.call(cbind, lapply(input, function(x) x[common,]))
}

#' subset beta value matrix by probes with minimal
#' number of non-NA.
#' 
#' @param betas beta value matrix
#' @param min_nonna_frac minimum fraction of Non-NA
#' @return subsetted beta value matrix
#' @examples
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas <- bSubAnyNonNAW(betas)
#' @export
bSubNonNAW <- function(betas, min_nonna_frac=0.5, max_na_cnt=NULL) {
    if (!is.null(max_na_cnt)) {
        return(betas[rowSums(is.na(betas)) <= max_na_cnt, ])
    }
    if (is.null(dim(betas))) { # should also work for vector
        betas[sum(!is.na(betas)) > min_nonna_frac*length(betas)]
    } else {
        betas[rowSums(!is.na(betas)) > min_nonna_frac*ncol(betas),]
    }
}

bSubCpGW <- function(betas) {
    betas[grep('^cg', rownames(betas)),]
}

bSubNoMaskW <- function(betas, platform='EPIC', refversion='hg38') {
    platform = 'EPIC'
    refversion = 'hg38'
    mft <- sesameDataGet(sprintf('%s.%s.manifest', platform, refversion))
    p <- intersect(rownames(betas), names(mft[!mft$MASK_general]))
    betas[p,]
}

cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
    cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
        f_row, f_col))
    cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
    namtx = is.na(mtx)
    good_row = rowSums(namtx) <= ncol(mtx) * (1-f_row)
    good_col = colSums(namtx) <= nrow(mtx) * (1-f_col)
    cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
    mtx[good_row, good_col]
}

sampleNMax <- function(df, column, k) {
    do.call(rbind, lapply(split(df, df[[column]]), function(df1) {
        if(nrow(df1) > 0) {
            sample_n(df1, min(nrow(df1),k))
        } else {
            NULL
        }
    }))
}

## sample members in groups
sampleGroup <- function(df, column, groups, k) {
    idx = df[[column]] %in% groups
    df_t = sample_n(df[idx,], k)
    df_n = df[!idx,]
    rbind(df_t, df_n)
}

mapMatrixW <- function(mat, val2num) {
    matrix(val2num[mat], nrow=nrow(mat), ncol=ncol(mat),
        dimnames=list(rownames(mat), colnames(mat)))
}

AFtoGenotype <- function(AFs) {
    AFs = matrix(cut(AFs, breaks=c(0,0.3,0.7,1), labels=c("R","H","A"), include.lowest=T), nrow=nrow(AFs), ncol=ncol(AFs), dimnames=list(rownames(AFs), colnames(AFs)))
}

#' calculate the AUC for scores to distinguish labels
#' 
#' @param labels 0-1 vector for case and control
#' @param scores a numerical vector of the same length as labels
#' @return the AUC (normalized Mann Whitney U)
auc_wmw2 = function(labels, scores) {
    labels <- as.logical(labels)
    n1 <- sum(labels)
    n2 <- sum(!labels)
    R1 <- sum(rank(scores)[labels])
    U1 <- R1 - n1 * (n1 + 1)/2
    U1/(n1 * n2)
}

#' sequentially remove duplicate from a list of char vector
removeDup = function(chars_list) {
    existing = NULL
    lapply(chars_list, function(chars) {
        chars = chars[!(chars %in% existing)]
        existing <<- c(existing, chars)
        chars
    })
}


imputeRowMean = function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    mtx
}

annotateW = function(accession, probeIDs, platform="MM285") {
    anno = tbk_data(sprintf("~/references/InfiniumArray/MM285/annotations/%s.tbk", accession), probes=probeIDs)
    anno[probeIDs,]
}


bothClusterSE = function(se, nrow_max = 3000, ncol_max = 3000) {
    if (nrow(se) > nrow_max || ncol(se) > ncol_max) {
        stop(sprintf("Too many rows (%d, max: %d) or columns (%d, max: %d). Abort.",
            nrow(se), nrow_max, ncol(se), ncol_max))
    }
    ord = both.cluster(assay(se))
    se[ord$row.clust$order, ord$column.clust$order]
}

buildSE = function(betas, probedf, meta, tissue_color, branch_color) {

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

meta2SE = function(meta, tbk_paths=NULL) {
    if (length(tbk_paths) == 0) {
        tbk_paths = NULL # TODO work this out
        ## tbk_paths ought to follow the same order as meta
    }

    betas = tbk_data(tbk_paths, max_pval=0.2)
    SummarizedExperiment(assays=list(betas=betas), colData=meta)
}

visualizeTissueSE = function(se, color=c("blueYellow","jet")) {
    ## this is the better version of sesame:::reference_plot_se(NULL, se)
    color = match.arg(color)
    if (color == "blueYellow") stop.points = c("blue","yellow")
    else stop.points = NULL
    cd = as_tibble(colData(se))
    rd = as_tibble(rowData(se))
    md = metadata(se)
    g = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points,
        dmin=0, dmax=1), name="b1", xticklabels=T, xticklabels.n=ncol(se))
    ## branch color bar (vertical)
    g = g + WColorBarV(rd$branch, RightOf("b1", width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh",
        yticklabels=T, yticklabel.side='l',
        label.space=0.01,
        label.use.data=TRUE, label.pad=1.5)
    ## tissue color bar (horizontal)
    g = g + WColorBarH(cd$branch, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$branch_color), name="ti",
        xticklabels=T, xticklabel.side='t',
        ## xticklabel.space=0.036,
        label.space=0.01,
        label.use.data=TRUE, label.pad=1.5)
    ## ## legends
    ## g = g + WLegendV("ti", TopRightOf("bh", just=c('left','top'), h.pad=0.5),
    ##     height=0.05)
    ## g = g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06, mar.top=0.1)
}

branchInferProbes = function(betas, meta, branch,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    branch_grouping = meta[[branch]]

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
    dfHypo = data.frame(delta_beta = delta_beta, auc = auc[names(delta_beta)],
        in_na = in_na[names(delta_beta)], out_na = out_na[names(delta_beta)],
        Probe_ID = names(delta_beta), branch=branch, type="Hypo")
    
    message(sprintf("Processed branch %s", branch))
    branchScanDeltaHypo(dfHypo)
    branchScanDeltaHype(dfHype)
    rbind(dfHypo, dfHype)
}

branchScanDeltaHypo = function(df, auc_min = 0.99) {
    message(sprintf("Found %d NA in Hypo", sum(is.na(df$delta_beta))))
    for (i in seq(-0.9,-0.1,by=0.05)) {
        probes = df$Probe_ID[(!is.na(df$delta_beta)) & df$delta_beta <= i & df$auc >= auc_min]
        message(sprintf("delta: %f Found %d probes", i, length(probes)))
    }
}

branchInferProbesMulti = function(betas, meta,
    hyper=FALSE, scan_delta = FALSE, scan_auc = FALSE) {

    meta$branch
    branch2samples = split(seq_along(meta$branch), meta$branch)
    ubs = do.call(cbind, lapply(branch2samples, function(samples) rowMaxs(betas[,samples])))
    lbs = do.call(cbind, lapply(branch2samples, function(samples) rowMins(betas[,samples])))

    group = do.call(rbind, lapply(seq_len(nrow(betas)), function(i) {
        ubs1 = ubs[i,]
        lbs1 = lbs[i,]
        if (all(is.na(ubs1)) || all(is.na(lbs1))) {
            hlb = NA; lub = NA;
        } else {
            hlb = max(lbs1, na.rm = T) # high group's lb
            lub = min(ubs1, na.rm = T) # low group's ub
        }
        is_l = ubs1 < hlb - 0.3
        is_h = lbs1 > lub + 0.3
        c(delta_beta = hlb - lub,
            ifelse(is.na(is_l) | is.na(is_h) | (is_l & is_h), 2,
                ifelse(is_l, 0, ifelse(is_h, 1, 2))))
    }))
    rownames(group) = rownames(betas)
    SummarizedExperiment(
        assays=list(membership=group[,2:ncol(group)]),
        rowData=data.frame(Probe_ID=rownames(group), delta_beta = group[,1]))
}

branchScanDeltaHype = function(df, auc_max = 0.01) {
    message(sprintf("Found %d NA in Hyper", sum(is.na(df$delta_beta))))
    for (i in seq(0.1,0.9,by=0.05)) {
        probes = df$Probe_ID[(!is.na(df$delta_beta)) & df$delta_beta >= i & df$auc <= auc_max]
        message(sprintf("delta: %f Found %d probes", i, length(probes)))
    }
}

branchInferProbes2 = function(res, delta_max = -0.6, auc_min = 0.99) {
    probes = names(which(res$delta_betas <= delta_max & res$auc >= auc_min))
    message(sprintf("Found %d probes", length(probes)))
    probes
}

branchInferProbes2Top = function(res, auc_min = 0.99, n=50) {
    probes = names(head(sort(res$delta_betas[res$auc >= auc_min]), n=n))
    message(sprintf("Found %d probes (max delta: %f)", length(probes), max(res$delta_betas[probes])))
    probes
}

## saveBranchProbes = function(branch, probes, date="NA", id_dir="~/Dropbox/Documents/paper/Ongoing_epigenetic_reconstruction/20210421_Branches_Human_EPIC/") {
##     write.table(tibble(probeID=probes, branchID=branch, date=date), file=sprintf("%s/%s.tsv", id_dir, branch), row.names=FALSE, quote=FALSE, sep="\t")
## }

readBranchProbes = function(id_dir = "~/Dropbox/Documents/paper/Ongoing_epigenetic_reconstruction/20210421_Branches_Human_EPIC/") {
    if (length(list.files(id_dir, ".rds"))==0) return(list())
    branches = sub(".rds","",sub("fit_","",list.files(id_dir,".rds")))
    do.call(rbind, lapply(branches, function(branch) {
        a = readRDS(sprintf("%s/fit_%s.rds", id_dir, branch))
        data.frame(delta_beta = a$delta_beta, auc = a[["auc"]][names(a$delta_beta)], branch = branch) %>% dplyr::filter(delta_beta <= -0.3, auc > 0.99) %>% rownames_to_column("Probe_ID") %>% dplyr::filter(grepl("^cg", Probe_ID))
    }))
}

## readBranchProbes0 = function(id_dir = "~/Dropbox/Documents/paper/Ongoing_epigenetic_reconstruction/20210421_Branches_Human_EPIC/", auc_min = 0.99, delta_max=-0.6, top_n=50) {
##     if (length(list.files(id_dir, ".rds"))==0) return(list())
##     fit_fns = list.files(id_dir, ".rds")
##     branchIDs = sub(".rds","",basename(fit_fns))
##     branchIDs = sub("fit_","",branchIDs)
##     setNames(lapply(fit_fns,
##         function(x) {
##             res = readRDS(sprintf("%s/%s", id_dir, x))
##             if (is.null(top_n)) {
##                 names(which(res$delta_betas <= delta_max & res$auc >= auc_min))
##             } else {
##                 names(head(sort(res$delta_betas[res$auc >= auc_min]), n=n))
##             }
##         }
##     ), branchIDs)
## }

## before regression check whether you have all the levels
checkLevels = function(data, meta) {
    apply(data, 1, function(dt) {
        all(sapply(split(dt, meta), function(x) sum(!is.na(x))>0))
    })
}

kcLoad = function(accession, home = "~/server/kyCG") {
    readRDS(sprintf("%s/%s.rds", home, accession))
}

#' @param priority database priority
#' @param home knowYourCG home
#' @return meta data tibble
#' @export
kcMeta = function(priority = 1, home = "~/server/kyCG", release = 1) {
    meta = read_csv(sprintf("%s/RELEASE_%d.csv", home, release))
    meta[meta$Priority <= priority,]
}

#' annotate CpGs using databases
#'
#' @param probes CGs to be annotated
#' @param accession database accession
#' @return character vector
#' @examples
#' anno = kcAnnotate(probes)
#' @export
kcAnnotate = function(probes,
    accession = "20210210_MM285_chromHMMconsensus",
    home = "~/server/kyCG",
    sep = ";") {
    
    db = kcLoad(accession, home)
    stopifnot(is(db, "list"))
    stopifnot(all(vapply(db, function(x) is(x, "character"), logical(1))))
    apply(do.call(cbind, lapply(names(db), function(nm) {
        ifelse(probes %in% db[[nm]], nm, NA); })), 1,
        function(x) paste(na.omit(x), collapse=sep))
}

#' calculate the mean beta value
#'
#' 
#' @param betas the beta value vector or matrix
#' @param accessions database accession, can be multiple
#' @return a tibble with mean betas
#' @export
kcMean = function(betas,
    accessions = "20210210_MM285_chromHMMconsensus",
    home = "~/server/kyCG") {
    
    do.call(rbind, lapply(accessions, function(acc) { # multiple accessions
        db = kcLoad(acc, home)
        do.call(rbind, lapply(names(db), function(nm) {
            probes = db[[nm]]
            if (is(betas, "matrix")) {
                tibble(
                    betaMean = colMeans(betas[probes, ,drop=FALSE], na.rm=TRUE),
                    sample = colnames(betas), db = nm)
            } else {
                tibble(betaMean = mean(betas[probes], na.rm=TRUE), db = nm)
            }
        }))
    }))
}
