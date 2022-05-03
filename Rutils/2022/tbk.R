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
