# SaveKaas = function(IDstring) { #save or update current kaas to schuifkaas, under name IDstring
#   if (!is.character(IDstring))
#     stop("a set k's are stored under a unique name which you should provide as parameter")
#   if (IDstring == "Current")
#     stop("Current is a protected word in this context, meaning the current public$kaas")
#   private$SchuifKaas[[IDstring]] <- private$SBkaas
# 
#' @name DiffKaas
#' @description extra return a diff of two sets of k's
#' @param ThisKaas
#' @param OtherKaas 
#' @param withProcess include its name
#' global, not to disrupt a logical object structure
#' @export
DiffKaas <- function(ThisKaas, OtherKaas, withProcess = T){ 
  #kaases should be dataframe-like with kaas properties, like SBcore$kaas 
  stopifnot("data.frame" %in% class(ThisKaas))
  stopifnot("data.frame" %in% class(OtherKaas))
  #i and j (indices of the boxes) might not consistent over worlds; remove
  kaasproperties <- c("i", "j", "k", "process", "fromSpecies", "fromSubCompart",
                      "fromScale", "fromAbbr", "toSpecies", "toSubCompart", "toScale", "toAbbr")
  stopifnot(all(names(ThisKaas) %in% kaasproperties))
  stopifnot(all(names(OtherKaas) %in% kaasproperties))
  ThisKaas$i <- NULL
  ThisKaas$j <- NULL
  OtherKaas$i <- NULL
  OtherKaas$j <- NULL
  
  #add distictive names
  ThisKaas$Origine <- "This"
  OtherKaas$Origine <- "Other"
  VolleKaas <- rbind(ThisKaas, OtherKaas)
  if (withProcess) {
    idVar <- c("process") #?
  } else {
    #remove process from dataframe, sum their kaas (each state-transition can have multiple processes)
    RHSformula <- c(kaasproperties[!kaasproperties %in% c("process", "kaas")],
                    "Origine") %>% 
      paste(sep = "+")
    AggFormula = as.formula(paste("kaas ~", RHSformula))
    VolleKaas <- aggregate(AggFormula, data = VolleKaas, FUN = sum)
    #? idVar <- c("i","j")
  }
  Wider <- pivot_wider(data = VolleKaas, values_from = c("k"),
                   names_from = c("Origine"))
}
