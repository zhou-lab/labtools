ggVerticalX <- function() {
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

ggPercentage <- function() {
    scale_y_continuous(labels=scales::percent)
}

ggNoLegend <- function() {
    guides(fill="none")
    ## or theme(legend.position="none")
    ## or scale_fill_discrete(guide="none")
}

view <- function(data, autofilter=TRUE) {
    ## source: https://jeromyanglim.tumblr.com/post/33825729070/function-to-view-r-data-frame-in-spreadsheet
    ## data: data frame
    ## autofilter: whether to apply a filter to make sorting and filtering easier
    open_command <- switch(Sys.info()[['sysname']],
        Windows= 'open',
        Linux  = 'xdg-open',
        Darwin = 'open')
    require(XLConnect)
    temp_file <- paste0(tempfile(), '.xlsx')
    wb <- loadWorkbook(temp_file, create = TRUE)
    createSheet(wb, name = "temp")
    writeWorksheet(wb, data, sheet = "temp", startRow = 1, startCol = 1)
    if (autofilter) setAutoFilter(wb, 'temp', aref('A1', dim(data)))
    saveWorkbook(wb, )
    system(paste(open_command, temp_file))
}

xtab_set <- function(A,B,AnotB = FALSE, BnotA = FALSE, AandB = FALSE){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B

    if (AnotB) {
        both[!(both %in% B)]
    } else if (BnotA) {
        both[!(both %in% A)]
    } else if (AandB) {
        both[(both %in% A) & (both %in% B)]
    } else {
        table(inA,inB)
    }
}
