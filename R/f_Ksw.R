#'@title soil water partitioning coefficient
#'@name f_Ksw
#'@description assess the Ksw if it's not given
#'@param ChemClass Class of chemical, see
#'@param Kow octanol water partitioning coefficient
#'@param CorgStandard Standard mass FRACTION organic carbon in soil/sediment [-]
#'@param a see QSAR table 
#'@param b see QSAR table 
#'@param RHOsolid density of "solids"
#'@param alt_form for base situations
#'@param Ksw_orig soil water partitioning coefficient as present in the data
#'@return Ksw
#'@export
f_Ksw <- function(Kow, pKa, CorgStandard , a, b, ChemClass, RHOsolid, alt_form, Ksw_orig){

  ifelse(alt_form,
         # TRUE, so the alt_form
         switch(ChemClass,
                "acid" = 10^(0.11*log10(Kow)+1.54) * CorgStandard * RHOsolid / 1000,
                "base" = 10^(pKa^0.65*(Kow/(1+Kow))^0.14) * CorgStandard * RHOsolid / 1000,
                #else
                {a * Kow^b * CorgStandard * RHOsolid / 1000}      
         ),
         # FALSE, NB not the alt_form
      switch(ChemClass,
          "acid" = 10^(0.54*log10(Kow)+1.11) * CorgStandard * RHOsolid / 1000,
          "base" = 10^(0.37*log10(Kow)+1.7) * CorgStandard * RHOsolid / 1000,
          "metal" = stop("Ksw Should be in the data"),
          "particle" = stop("Ksw Should be in the data"),
          #else
          {a * Kow^b * CorgStandard * RHOsolid / 1000}
      )
    )
}
