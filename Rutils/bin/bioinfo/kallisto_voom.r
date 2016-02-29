#!/usr/bin/env r
## kallisto_voom.r -g mm10 -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o 
suppressMessages(library(docopt))

"Usage:
  kallisto_voom.r -a CONDITION1 -b CONDITION2 -A SNAMES1 -B SNAMES2 -o OUTPUT

Options:
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A SNAMES1      sample names for condition1 (comma-separated, kallisto/sname/abundance.tsv)
  -B SNAMES2      sample names for condition2 (comma-separated, kallisto/sname/abundance.tsv)
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
sink(stderr())
# cat(str(opt))

cat(format(Sys.time()), "Read samples...\n")
snames1 <- unlist(strsplit(opt$A,','))
snames2 <- unlist(strsplit(opt$B,','))
n1 <- length(snames1)
n2 <- length(snames2)
sample_list_n <- c(snames1, snames2)
for (i in 1:length(sample_list_n)) {
    tmp = read.table(file=paste0("kallisto/",sample_list_n[i],"/abundance.tsv"), header=T)
    assign(sample_list_n[i], tmp)       # attach each data frame to the current environment
}

sample_list = mget(sample_list_n)

## give the list unique names 
sample_list_uni = Map(function(x, i) setNames(x, ifelse(names(x) %in% "target_id",
    names(x), sprintf('%s.%d', names(x), i))), sample_list, seq_along(sample_list))

full_kalli = Reduce(function(...) merge(..., by = "target_id", all=T), sample_list_uni)
 
tpm_vals = full_kalli[, grep("tpm", names(full_kalli))]
rownames(tpm_vals) = full_kalli$target_id

## then make a contrast matrix
groups = rep(c(opt$a, opt$b), c(n1,n2))
condition <- model.matrix(~0 + groups)
colnames(condition) = c(opt$a, opt$b)
suppressPackageStartupMessages(library(limma))
cont_matrix = makeContrasts(paste0(opt$a,"-",opt$b), levels = condition)

## this spits out an EList object 
v = voom(counts = tpm_vals, design = condition)
fit = lmFit(v, condition)
fit = contrasts.fit(fit, cont_matrix)
fit = eBayes(fit)
top_table = topTable(fit, n = 100000000, sort.by = "p")

colnames(tpm_vals) <- names(sample_list)

write.table(merge(tpm_vals,top_table,by=0), file=opt$o, sep="\t", quote=F, row.names=F)
cat(format(Sys.time()), "Done.\n")
