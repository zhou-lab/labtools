
getInfiniumProbesByCoord <- function(Probe_ID2, platform="HM450") {
    mft = sesameData_getManifestGRanges(platform)
    data.frame(Probe_ID2=Probe_ID2,Probe_ID=names(mft)[match(Probe_ID2,paste0(seqnames(mft),"_",start(mft)))])
}

convertToArray <- function(se, platform="HM450") {
    mft = sesameData_getManifestGRanges(platform, genome="hg38")
    idx = match(rowData(se)$Probe_ID, paste0(seqnames(mft),"_",start(mft)))
    se1 = se[!is.na(idx),]
    rowData(se1)$Probe_ID = names(mft)[na.omit(idx)]
    se1
}

cgNumberToCoord <- function(cgNumbers, platform=NULL) {
    ## cgNumberToCoord(c("cg14640477", "cg01384075", "cg15835644"))
    ## [1] "chr10_26457265"  "chr10_102002633" "chr7_26360171"
    ## chr10_26457265, chr10_102002633, chr7_26360171
    if (is.null(platform)) {
        platform = inferPlatformFromProbeIDs(cgNumbers)
    }
    mft = sesameData_getManifestGRanges(platform)
    cgNo0 = names(mft)
    coords = paste0(seqnames(mft),"_",start(mft))
    coords[match(cgNumbers, cgNo0)]
}
