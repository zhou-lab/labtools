
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

rowMax <- function(x) {apply(x,1,max)}

suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(devtools)))
## suppressWarnings(suppressPackageStartupMessages(library(wheatmap)))
## suppressWarnings(suppressPackageStartupMessages(library(sesame)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr, quiet=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges, quiet=TRUE)))
## suppressWarnings(suppressPackageStartupMessages(library(limma)))
## theme_wz <- theme_set(theme_classic(15))
## theme_wz <- theme_update(
##   axis.line.x = element_line(colour = "black"),
##   axis.line.y = element_line(colour = "black"),
##   axis.text = element_text(colour='black'))

theme_wz_light <- theme_set(theme_linedraw(15))
theme_wz_light <- theme_update(
  ## text = element_text(face='bold'),
  ## axis.line.x = element_line(colour = "black", size=1.2),
  ## axis.line.y = element_line(colour = "black", size=1.2),
  axis.text = element_text(colour='black', size=15),
  panel.border = element_rect(linetype = "solid", colour = "black", size=1.2),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank())

theme_wz <- theme_set(theme_linedraw(16))
theme_wz <- theme_update(
  text = element_text(face='bold'),
  ## axis.line.x = element_line(colour = "black", size=1.2),
  ## axis.line.y = element_line(colour = "black", size=1.2),
  axis.text = element_text(colour='black', size=16),
  panel.border = element_rect(linetype = "solid", colour = "black", size=1.2),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
