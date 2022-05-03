readExcelColors <- function(path) {
    sheets = excel_sheets(path)
    sheets = sheets[endsWith(sheets, "colors")]
    colors = lapply(sheets, function(sheet) read_excel(path, sheet) %>% with(setNames(HEX, LABEL)))
    names(colors) = sheets
    colors
}
