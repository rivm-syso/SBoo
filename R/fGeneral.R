# fGeneral

#' @description 2 basic ordinary differential equation (ODE) functions for simplebox; the function for the ode-call;
#' see desolve/rootsolve packages
#' @param t time (vector ?)
#' @param m  (i) = initial mass
#' @param parms = (K, e) i.e. the matrix of speed constants and the emissions as vector, or functions
#' @returns dm (i) = change in mass as list

SimpleBoxODE = function(t, m, parms) {
  dm <- with(parms, 
             K %*% m + e)
  return(list(dm, signal = parms$e)) 
}

SimpleBoxODEapprox = function(t, m, parms) {
  # define what parms are needed
  with(as.list(c(parms, m)), {
    e <- c(rep(0, length(SBNames)))
    
    for (name in names(parms$emislist)) {
      e[grep(name, SBNames)] <- parms$emislist[[name]](t) 
    }
    dm <- with(parms, K %*% m + e) 
    return(list(dm, signal = e))
  })
}

#' @description Steady state solver function
#' @param k the first order rate constant matrix
#' @param m emission vector
#' @param parms = empty list (needed so that the SteadyStateSolver and DynamicSolver can be called in the same manner from the SolverModule)
#' @returns dm (i) = change in mass as list

SteadyStateSolver <- function(k, e, parms){
  SBNames = colnames(k)
  tmax=1e20 # solution for >1e12 year horizon
  dm <- rootSolve::runsteady(
    y = rep(0,nrow(k)),
    times = c(0,tmax),
    func = SimpleBoxODE,
    parms = list(K = k, 
                 e = e)
  )
  return(dm)
}

#' @description Dynamic solver function
#' @param k the first order rate constant matrix
#' @param m list of emission approx functions per compartment
#' @param parms = list containing nTIMES, tmin and tmax
#' @returns a list containing main (a matrix with the masses per compartment per
#'  timestep) and signals (a matrix with the emissions per compartment per timestep)

DynamicSolver <- function(k, e, parms) {
    if(is.null(parms$rtol_ode)){
    rtol_ode = 1e-11
  } else rtol_ode = parms$rtol_ode
  if(is.null(parms$atol_ode)){
    atol_ode = 1e-3
  } else atol_ode = parms$atol_ode
  tmax = parms$tmax 
  nTIMES = parms$nTIMES 
  SB.K = k
  SBNames = colnames(k)
  SB.m0 = rep(0, length(SBNames))
  tmin = parms$tmin
  SBtime <- seq(tmin,tmax,length.out = nTIMES)
  
  out <- deSolve::ode(
    y = as.numeric(SB.m0),
    times = SBtime,
    func = SimpleBoxODEapprox,
    parms = list(K = SB.K, 
                 SBNames=SBNames,
                 emislist = e),
   rtol = rtol_ode, atol =  atol_ode)

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

#' @description helper function for Make_inv_unif01
triangular_cdf_inv = function(u, # LH scaling factor
                              a, # Minimum
                              b, # Maximum
                              c) { # Peak value
  ifelse(u < (c-a)/(b-a),
         a + sqrt(u * (b-a) * (c-a)),
         b - sqrt((1-u) * (b-a) * (b-c)))
}

#create a function for transformation of lhs range (0-1) to actual variable range (inverse of the 0-1 cdf)
Make_inv_unif01 = function(fun_type = "triangular", pars) {
  if (!fun_type %in% c("triangular", "normal", "uniform")) {
    stop("! fun_type %in% c('triangular', 'normal', 'uniform')")
  }
  if (fun_type == "triangular") {
    if (!(inherits(pars, "list") && length(pars) == 3)) {
      stop(
        "the triangular is created using a list of three parameters, a = minimum, b = maximum, c = peak")
    }
    a <- pars[["a"]]
    b <- pars[["b"]]
    c <- pars[["c"]]
    return(function(x) {
      triangular_cdf_inv(x, a, b, c)
    })
  }
  if (fun_type == "normal") {
    if (!(inherits(pars, "list")) && length(pars) == 2) {
      stop("the normal is created using a list of two parameters, a = mean, b = sigma, c = peak")
    }
    mu <- pars[["a"]]
    sig <- pars[["b"]]
    return(function(x) {
      qnorm(p = x, mean = mu, sd = sig)
    })
  }
  if (fun_type == "uniform") {
    if (!(inherits(pars, "list")) && length(pars) == 2) {
      stop("the uniform is created using a list of two parameters, a = minimum, b = maximum")
    }
    minx <- pars[["a"]]
    maxx <- pars[["b"]]
    return(function(x) {
      minx + x * (maxx - minx)
    })
  }
  
}