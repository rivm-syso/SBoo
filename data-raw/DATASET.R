## code to prepare `DATASET` dataset goes here

The3D <- c("Scale","SubCompart","Species")
usethis::use_data(The3D, internal = FALSE, overwrite = T)
#initial dataframes from M 

Defs = c(
  "ScaleSubCompartData", #Bevat landFRAC
  "ScaleSpeciesData",
  "SubCompartSpeciesData",
  "ScaleSheet",
  "SubCompartSheet",
  "SpeciesSheet",
  "ScaleProcesses",
  "SubCompartProcesses",
  "SpeciesProcesses",
  "Compartments",
  "Substances",
  "SpeciesCompartments",
  "SubstanceCompartments",
  "SubstanceSubCompartSpeciesData",
  "CONSTANTS",
  "MatrixSheet",
  "FlowIO",
  "QSARtable",
  "SomeFromTo",
  "Units"
)
usethis::use_data(Defs, internal = FALSE, overwrite = T)

RowIdentifyers = c(
  The3D, "Species", 
  "to.Scale", "to.SubCompart", "to.Species", 
  "Substance", "process", "Matrix", "from", "to", "QSAR.ChemClass"
)

usethis::use_data(RowIdentifyers, internal = FALSE, overwrite = T)