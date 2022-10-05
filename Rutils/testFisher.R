args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
options(scipen=99)
options(digits=3)
df = read.table(args[1], header=T) %>% dplyr::mutate(
    nf=nfmq+nfq,
    nqmf=nq-nfq,
    numfq=nu-nf-nq+nfq,
    dep = phyper(nfq, nf, nu-nf, nq, lower.tail=TRUE),
    enrich = phyper(nfq-1, nf, nu-nf, nq, lower.tail=TRUE),
    odds_ratio=(nfq/nfmq)/(nqmf/numfq))
write.table(format(df), file=args[2], quote=FALSE, row.names=FALSE, sep="\t")
