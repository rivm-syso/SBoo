#' @title ConcentrationModule
#' @description module for calculation of concentrations from masses
#' @param TheCore centre "World" object
#' @param input any other input that can be given, can be used to force the function to use a dataframe
#' @import R6
#' @export
ConcentrationModule <- 
  R6::R6Class("ConcentrationModule",
      public = list(
        initialize = function(TheCore, input){
          #browser()
          private$MyCore <- TheCore
          
        },
        #'@description Leading function to obtain concentration from other world parameters
        #'@input all input is taken from the core
        GetConc = function(input){
          #browser()
          solution <- private$MyCore$Solution()
          states <- private$MyCore$states$asDataFrame
          volume <- private$MyCore$fetchData("Volume")
          fracw <- private$MyCore$fetchData("FRACw")
          fraca <- private$MyCore$fetchData("FRACa")
          rho <- private$MyCore$fetchData("rhoMatrix")
          
          longsolution <- solution |>
            left_join(states, by = "Abbr") |>
            left_join(private$MyCore$fetchData("Volume"), by = c("SubCompart", "Scale")) |>
            mutate(Concentration = EqMass / Volume)
          conctobecor <- longsolution  |>
            filter(SubCompart %in% c("marinesediment", "freshwatersediment", "lakesediment", 
                                     "agriculturalsoil", "naturalsoil", "othersoil"))
          conctobecor <- conctobecor |>
            left_join(fracw, by = c("SubCompart", "Scale")) |>
            left_join(fraca, by = c("SubCompart", "Scale")) |> 
            left_join(rho, by = c("SubCompart")) |> 
            mutate(rhowat = 998)
          
          conctobecor <- conctobecor |> 
            mutate(
              AdjustedConc = private$adjustconc(Concentration, FRACw, FRACa, rhoMatrix, rhowat)
            )
          
          longsolution <- longsolution |>
            left_join(
              conctobecor |> select(Scale, SubCompart, Species, AdjustedConc), 
              by = c("Scale", "SubCompart", "Species")
            ) |>
            mutate(
              Concentration = coalesce(AdjustedConc, Concentration)
            ) |>
            select(-AdjustedConc)
          
          FinalConcentration <- longsolution |>
            select(Abbr, Concentration)
        
          return(FinalConcentration)
        
        }
      ),
        # concentrationcalc = function() 
      private = list(
        MyCore = NULL, 
        #'@description helper function to adjust concentration for soil and sediment 
        #'@param CompConc Concentration in each compartment [mass/volume]
        #'@param FRACw the fraction of water in the matrix
        #'@param FRACa the fraction of air in the matrix 
        #'@param RHOsolid the density of the solid phase of the material [kg m-3]
        #'@param RhoWater_value the density of water [kg m-3]
        adjustconc = function(CompConc, Fracw, Fraca, RHOsolid, RhoWater_value){
            CompConc  / (Fracw * RhoWater_value + (1 - Fracw - Fraca) * RHOsolid)
          }
      )
  )

