# fGeneral

#' @name  getConst
#' @description grab from the web: expand ... data.frames
#' global, not to disrupt a logical object structure
#' @export
getConst <- function(symbol) {
  res <- constants::codata$value[constants::codata$symbol == symbol]
  if (length(res) == 0 || is.na(res)) stop (paste(symbol, "not found in constants::codata, use ConstGrep()"))
  return(res)
}

#' @name  ConstGrep
#' @description wrapper around constants::codata
#' @param grepSearch term to search for
#' @param ... passed on to grep, use for instance ConstGrep("Gas", ignore.case=TRUE)
#' @export
ConstGrep <- function(grepSearch, ...){
  constants::codata[grep(grepSearch, constants::codata$quantity, ...), ]
}

#' @name  expand.grid.df
#' @description grab from the web: expand ... data.frames
#' @param ... {n} data.frame
#' global, not to disrupt a logical object structure
#' @export
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

#' @name  dframe2excel
#' @description grab from the web: expand ... data.frames
#' @param dframe data to output to 
#' @param outxlsx filename
#' global, not to disrupt a logical object structure
#' @export
dframe2excel <- function(dframe, outxlsx = "sb2excel.xlsx") {
  if(!endsWith(outxlsx, ".xlsx")) {
    outxlsx <- paste0(outxlsx, ".xlsx")
  }
  if (file.exists(outxlsx)) {
    existSheets <- openxlsx::getSheetNames(file = outxlsx)
    wb <- openxlsx::loadWorkbook(file = outxlsx)
  } else {
    wb <- openxlsx::createWorkbook()
    existSheets <- c("")
  }
  sheetName <- attr(dframe, which = "sheetName")
  if (length((sheetName)==0) | !is.character(sheetName)) {
    sheetName <- do.call(paste0,as.list(names(dframe)))
    if (nchar(sheetName) > 28) {
      #maximaze length to 28 chrs
      prtNames <- names(dframe)[1:min(28,length(names(dframe)))]
      charsPName <- floor(28 / length(prtNames))
      dfNames <- lapply(names(dframe), function(x) {substr(x,1,charsPName)})
      sheetName <- substr(do.call(paste0,dfNames), start = 1, stop = 28)
    }
  } 
  if (sheetName %in% existSheets) 
    openxlsx::removeWorksheet(wb=wb, sheet = sheetName)
  openxlsx::addWorksheet(wb=wb, sheetName = sheetName)
  openxlsx::writeData(wb, sheet = sheetName, dframe)
  openxlsx::saveWorkbook(wb, outxlsx, overwrite = TRUE)
}



