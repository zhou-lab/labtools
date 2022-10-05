args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
options(scipen=99)
options(digits=3)
df = read.table(args[1], header=T) %>% dplyr::mutate(
    nf=nfmq+nfq,
    nqmf=nq-nfq,
    numfq=nu-nf-nq+nfq,
    p_val_enr = phyper(nfq, nf, nu-nf, nq, lower.tail=TRUE),
    p_val_dep = phyper(nfq-1, nf, nu-nf, nq, lower.tail=FALSE),
    odds_ratio=(nfq/nfmq)/(nqmf/numfq)) %>% dplyr::arrange(p_val_enr)
write.table(format(df), file=args[2], quote=FALSE, row.names=FALSE, sep="\t")
