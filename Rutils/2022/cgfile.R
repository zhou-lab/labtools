load_cgfile_to_matrix <- function(cg_file_path, probes, min_cov=5) {
    ## Example input:
    ## cg_file_path = "/mnt/isilon/zhou_lab/projects/20230727_all_public_WGBS/hg38/2023_Loyfer.cg"
    ## probes = c("chr8_28338729", "chr1_244941270", "chr2_25563815")
    temp_file <- tempfile(tmpdir="~/tmp/")
    write.table(data.frame(a=probes), file=temp_file,
        row.names = FALSE, col.names = FALSE, quote=FALSE)

    snames = read.table(sprintf("%s.idx", cg_file_path))$V1
    mtx = bind_cols(mclapply(snames, function(sname) read.table(text=system(sprintf("yame subset %s %s | yame rowsub -R ~/references/hg38/KYCGKB_hg38/cpg_nocontig.cr -L ~/tmp//file8d8811d97681b - | yame unpack -f %d -", cg_file_path, sname, min_cov), intern=TRUE), header=FALSE, sep="\t"), mc.cores=15))
    rownames(mtx) = probes
    colnames(mtx) = snames
    mtx
}
