

txdb <- makeTxDbFromGFF("/primary/vari/genomicdata/genomes/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz",format="gtf",circ_seqs=character())

ebg <- exonsBy(txdb, by="gene")

se <- summarizeOverlaps(features=ebg, reads=c("bam/YS_female_mut.bam", "bam/YS_female_wt.bam"),mode="Union",singleEnd=F, ignore.strand=F, fragments=T)


seqlevelsStyle(ebg) <- "UCSC"
fpkm()
