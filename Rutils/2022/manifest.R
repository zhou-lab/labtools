#' Annotate a data.frame using manifest
#'
#' @param df input data frame with Probe_ID as a column
#' @param probe_id the Probe_ID column name, default to "Probe_ID" or
#' rownames
#' @param platform which array platform, guess from probe ID if not given
#' @param genome the genome build, use default if not given
#' @return a new data.frame with manifest attached
#' @examples
#' df <- data.frame(Probe_ID = c("cg00101675_BC21", "cg00116289_BC21"))
#' attachManifest(df)
#' @export
attachManifest <- function(
    df, probe_id="Probe_ID", platform=NULL, genome=NULL) {
    stopifnot(is(df, "data.frame"))
    stopifnot(probe_id %in% colnames(df))

    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(df[[probe_id]]) }

    genome <- sesameData_check_genome(genome, platform)

    mft <- sesameDataGet(sprintf("%s.%s.manifest", platform, genome))
    cbind(df, as.data.frame(mft)[match(df[[probe_id]], names(mft)),])
}
