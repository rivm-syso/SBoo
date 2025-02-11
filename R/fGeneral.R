# fGeneral

#' @description 4 basic ODE functions for simplebox; the function for the ode-call;
#' see desolve/rootsolve packages
#' @param t time (vector ?)
#' @param m  (i) = initial mass
#' @param parms = (K, e) i.e. the matrix of speed constants and the emissions as vector, or functions
#' @returns dm (i) = change in mass as list
#' 
# SteadyODE <- function(k, m, tol=1e-30, parms) {
#   dm <- solve(k, -m, tol = tol)
#   return(list(dm))
# }

SimpleBoxODE = function(t, m, parms) {
  dm <- with(parms, K %*% m + e)
  return(list(dm, signal = parms$e)) 
}

SteadyODE <- function(k, m, parms){
  tmax=1e20 # solution for >1e12 year horizon
  sol <- rootSolve::runsteady(
    y = rep(0,nrow(k)),
    times = c(0,tmax),
    func = SimpleBoxODE,
    parms = list(K = k, e = m)
  )
}

EmisBoxODE <- function(t, m, parms) {
  e_t <- parms$e(t)
  dm <- with(parms, K %*% m + e_t)
  return(dm)
}

EventODE <- function(k, m, parms) {
  with(as.list(c(parms, m)), {
    dm <- k %*% m
  })
  return(list(dm))
}

# ODEapprox <- function(times, y, parms) {
#   
#   t <- times
#   with(as.list(c(parms, y)), {
#     e <- c(rep(0, length(SBNames)))
#     for (name in names(parms$emislist)) {
#       e[grep(name, SBNames)] <- parms$emislist[[name]](t)
#     }
#     dm <- with(parms, K %*% y + e) 
# 
#     return(list(dm, signal = e))
#   })
# }

ODEapprox = function(t, m, parms) {
  with(as.list(c(parms, m)), {
    e <- numeric(length(SBNames))  # Initialize emissions vector with zeros
    #print(e)
    for (name in names(parms$emislist)) {
      idx <- match(name, SBNames)  # Find the exact index of the compartment
      if (!is.na(idx)) {
        e[idx] <- parms$emislist[[name]](t)  # Assign emission to the correct compartment
      }
    }
    # Debugging: Print emissions vector
    #print(paste("Time:", t, "Emissions:", paste(e, collapse = ", ")))
    
    dm <- with(parms, K %*% m + e) 
    return(list(dm, signal = e))
  })
}

ApproxODE <- function(k, m, parms) {
  SBNames <- colnames(k)  # Assuming k is a matrix
  SB.m0 <- rep(0, length(SBNames))
  SBtime <- seq(0, parms$tmax, length.out = parms$nTIMES)
  
  out <- deSolve::ode(
    y = as.numeric(SB.m0),
    times = SBtime,
    func = ODEapprox,
    parms = list(K = k, SBNames = SBNames, emislist = m),
    rtol = 1e-30, atol = 1e-3
  )
  
  colnames(out)[1:length(SBNames) + 1] <- SBNames
  colnames(out)[grep("signal", colnames(out))] <- paste("signal", SBNames, sep = "2")
  
  signal_cols <- grep("^signal", colnames(out))
  
  # Extract the "signal" columns into a new matrix
  signal_matrix <- out[, signal_cols, drop = FALSE]
  
  # Remove the "signal" columns from the original matrix
  cols_to_exclude <- c(1, signal_cols)
  out <- out[, -cols_to_exclude, drop = FALSE]
  
  # Remove "signal" from the column names in the signal_matrix
  colnames(signal_matrix) <- sub("^signal2", "", colnames(signal_matrix))
  
  return(list(main = out, signals = signal_matrix))
}

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


#' Title
#' @param df input object df, or df-alike
#' @param callingName name of the calling function, or? user informative string
#' @param mustHavecols column names that the df should have
#' @return side effect (stop) only
#' @export
is.df.with <- function(df, callingName, mustHavecols, reallystop = T){
  if (! "data.frame" %in% class(df)) {
    rep_text <- paste("Not a data.frame-ish parameter in", callingName)
    if (reallystop){
      stop(rep_text)
    } else warning(rep_text)
  }
  missingCol <- mustHavecols[!mustHavecols %in% names(df)]
  if (length(missingCol) > 0) {
    rep_text <- do.call(paste, as.list(c("Missing column", missingCol, 
                                         "in", callingName)))
    if (reallystop){
      stop(rep_text)
    } else warning(rep_text)
  }
  
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
