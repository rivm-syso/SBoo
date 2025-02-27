#################### Functions for solution plots

# Solution plot for deterministic steady state output
DetSSPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the solution
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  if (!is.null(subcompart)) {
    
    solution <- solution[solution$SubCompart %in% subcompart, ]
  } else {
    solution <- solution
  }
  
  plot <- ggplot(solution, aes(x = SubCompart, y = Mass_kg, fill = SubCompart)) + 
    theme_bw() + 
    geom_col() +  # Use geom_col() to plot the actual values of Mass_kg
    labs(title = paste0("Steady state mass at ", scale, " scale"),
         x = "",
         y = paste0("Mass of substance [kg]")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) 
}

# Solution plot for deterministic dynamic output
DetDynSolPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the solution
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'time', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  if (!is.null(subcompart)) {
    solution <- solution[solution$SubCompart %in% subcompart, ]
  } else {
    solution <- solution
  }
  
  # Convert time from seconds to years
  solution$time <- as.numeric(solution$time)
  solution$Year <- solution$time / (365.25 * 24 * 3600)
  
  plot <- ggplot(solution, aes(x = Year, y = Mass_kg, group = SubCompart, color = SubCompart)) + 
    theme_bw() + 
    geom_line() +
    labs(title = paste0("Dynamic mass at ", scale, " scale"),
         x = "Year",
         y = paste0("Mass of substance [kg]")) +
    scale_y_continuous(labels = scales::label_scientific()) +
    guides(color = guide_legend(title = NULL))
}

# Solution plot for probabilistic dynamic output

ProbDynSolPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the solution
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'time', 'RUNs', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  if(length(subcompart) != 1){
    stop("Please select 1 subcompartment")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  if (!is.null(subcompart)) {
    solution <- solution[solution$SubCompart %in% subcompart, ]
  } else {
    solution <- solution
  }
  
  # Convert time from seconds to years
  solution$time <- as.numeric(solution$time)
  solution$Year <- solution$time / (365.25 * 24 * 3600)
  
  summary_stats <- solution |>
    group_by(Year) |>
    summarise(
      Mean_Value = mean(Mass_kg, na.rm = TRUE),
      SD_Value = sd(Mass_kg, na.rm = TRUE),
      Lower_CI = Mean_Value - 1.96 * SD_Value / sqrt(n()),
      Upper_CI = Mean_Value + 1.96 * SD_Value / sqrt(n())
    ) |>
    ungroup()
  
  ggplot(summary_stats, aes(x = Year, y = Mean_Value)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, fill = "blue") +
    labs(title = paste("Dynamic mean mass in", subcompart, " at ", scale, " scale"),
         subtitle = "with uncertainty bands over time",
         x = "Year",
         y = paste("Mass of substance in [kg]")) +
    theme_minimal()+
    guides(color = guide_legend(title = NULL))
}

# Solution plot for probabilistic steady state output
ProbSSSolPlot <- function(scale = NULL){
  
  # Get the solution
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'RUNs', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  ggplot(solution, aes(x = SubCompart, y = Mass_kg, fill = SubCompart)) +
    geom_violin()+
    theme_bw() +
    labs(title = paste("Steadt state mass at ", scale, "scale"),
         x = "",
         y = paste("Mass of substance [kg]")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_discrete() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

#################### Functions for concentration plots

# Concentration plot for deterministic dynamic output
DetSSConcPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the solution
  concentration <- merge(World$Concentration(), World$states$asDataFrame, by = "Abbr")
  concentration <- concentration[c('SubCompart', 'Scale', 'Species', 'Concentration', 'Unit')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(concentration$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(concentration$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(concentration)
  cnames <- cnames[!cnames %in% c("Species", "Concentration")]
  formula <- as.formula(paste("Concentration ~", paste(cnames, collapse = " + ")))
  concentration <- aggregate(formula, data = concentration, sum)
  
  if (!is.null(scale)) {
    concentration <- concentration[concentration$Scale %in% scale, ]
  } else {
    concentration <- concentration
  }
  
  if (!is.null(subcompart)) {
    concentration <- concentration[concentration$SubCompart %in% subcompart, ]
  } else {
    concentration <- concentration
  }
  
  concentration$SubCompartUnit <- paste0(concentration$SubCompart, " (", concentration$Unit, ")")
  
  plot <- ggplot(concentration, aes(x = SubCompartUnit, y = Concentration, fill = SubCompartUnit)) + 
    theme_bw() + 
    geom_col() +  
    labs(title = paste0("Steady state concentration at ", scale, " scale"),
         x = "",
         y = paste0("Concentration of substance")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) 
}

# Concentration plot for deterministic dynamic output
DetDynConcPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the solution
  concentration <- merge(World$Concentration(), World$states$asDataFrame, by = "Abbr")
  concentration <- concentration[c('SubCompart', 'Scale', 'Species', 'time', 'Concentration', 'Unit')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(concentration$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(concentration$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(concentration)
  cnames <- cnames[!cnames %in% c("Species", "Concentration")]
  formula <- as.formula(paste("Concentration ~", paste(cnames, collapse = " + ")))
  concentration <- aggregate(formula, data = concentration, sum)
  
  if (!is.null(scale)) {
    concentration <- concentration[concentration$Scale %in% scale, ]
  } else {
    concentration <- concentration
  }
  
  if (!is.null(subcompart)) {
    concentration <- concentration[concentration$SubCompart %in% subcompart, ]
  } else {
    concentration <- concentration
  }
  
  concentration$SubCompartUnit <- paste0(concentration$SubCompart, " (", concentration$Unit, ")")
  
  # Convert time from seconds to years
  concentration$time <- as.numeric(concentration$time)
  concentration$Year <- concentration$time / (365.25 * 24 * 3600)
  
  plot <- ggplot(concentration, aes(x = Year, y = Concentration, group = SubCompartUnit, color = SubCompartUnit)) + 
    theme_bw() + 
    geom_line() +
    labs(title = paste0("Dynamic concentration at ", scale, " scale"),
         x = "Year",
         y = "Concentration of substance") +
    scale_y_continuous(labels = scales::label_scientific()) +
    guides(color = guide_legend(title = NULL))
}

# Concentration plot for probabilistic dynamic output
ProbDynConcPlot <- function(scale = NULL, subcompart = NULL){
  
  # Get the concentration
  concentration <- merge(World$Concentration(), World$states$asDataFrame, by = "Abbr")
  concentration <- concentration[c('SubCompart', 'Scale', 'Species', 'time', 'RUNs', 'Concentration', 'Unit')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  if(length(subcompart) != 1){
    stop("Please select 1 subcompartment")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(concentration$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Make sure the selected subcompartments exist
  if (!is.null(subcompart) && !all(subcompart %in% unique(concentration$SubCompart))) {
    stop("One or more selected subcomparts do not exist")
  }
  
  # Aggregate over species
  cnames <- names(concentration)
  cnames <- cnames[!cnames %in% c("Species", "Concentration")]
  formula <- as.formula(paste("Concentration ~", paste(cnames, collapse = " + ")))
  concentration <- aggregate(formula, data = concentration, sum)
  
  if (!is.null(scale)) {
    concentration <- concentration[concentration$Scale %in% scale, ]
  } else {
    concentration <- concentration
  }
  
  if (!is.null(subcompart)) {
    concentration <- concentration[concentration$SubCompart %in% subcompart, ]
  } else {
    concentration <- concentration
  }
  
  unit <- unique(concentration$Unit)
  
  # Convert time from seconds to years
  concentration$time <- as.numeric(concentration$time)
  concentration$Year <- concentration$time / (365.25 * 24 * 3600)
  
  summary_stats <- concentration |>
    group_by(Year) |>
    summarise(
      Mean_Value = mean(Concentration, na.rm = TRUE),
      SD_Value = sd(Concentration, na.rm = TRUE),
      Lower_CI = Mean_Value - 1.96 * SD_Value / sqrt(n()),
      Upper_CI = Mean_Value + 1.96 * SD_Value / sqrt(n())
    ) |>
    ungroup()
  
  ggplot(summary_stats, aes(x = Year, y = Mean_Value)) +
    geom_line(color = "green", size = 1) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, fill = "green") +
    labs(title = paste0("Dynamic mean concentration in ", subcompart, " at ", scale, " scale"),
         subtitle = "with uncertainty bands over time",
         x = "Year",
         y = paste0("Concentration of substance [", unit, "]")) +
    theme_minimal()+
    guides(color = guide_legend(title = NULL))
}

# Concentration plot for probabilistic steady state output
ProbSSConcPlot <- function(scale = NULL){
  
  # Get the concentrations
  concentration <- merge(World$Concentration(), World$states$asDataFrame, by = "Abbr")
  concentration <- concentration[c('SubCompart', 'Scale', 'Species', 'RUNs', 'Concentration', 'Unit')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(concentration$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Aggregate over species
  cnames <- names(concentration)
  cnames <- cnames[!cnames %in% c("Species", "Concentration")]
  formula <- as.formula(paste("Concentration ~", paste(cnames, collapse = " + ")))
  concentration <- aggregate(formula, data = concentration, sum)
  
  if (!is.null(scale)) {
    concentration <- concentration[concentration$Scale %in% scale, ]
  } else {
    concentration <- concentration
  }
  
  concentration$SubCompartUnit <- paste0(concentration$SubCompart, " (", concentration$Unit, ")")
  
  ggplot(concentration, aes(x = SubCompartUnit, y = Concentration, fill = SubCompart)) +
    geom_violin()+
    theme_bw() +
    labs(title = paste0("Steady state concentration at ", scale, " scale"),
         x = "",
         y = paste0("Concentration")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_fill_discrete() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

######################## Mass distribution plot functions ######################

# For deterministic steady state masses
DetSSMassDist <- function(scale = NULL){
  
  # Get the solution and join to states df
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  ggplot(solution, aes(area = Mass_kg, fill = Mass_kg, label = SubCompart)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white",
                                  place = "centre",
                                  size = 15) +
    labs(title = paste0("Distribution of steady state masses at ", scale, " scale"),
         fill = "Mass [kg]") 
}

# For probabilistic steady state masses
ProbSSMassDist <- function(scale = NULL){
  
  # Get the solution and join to states df
  solution <- merge(World$Solution(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'RUNs', 'Mass_kg')]
  
  # Make sure 1 scale is selected
  if(length(scale) != 1){
    stop("Please select 1 scale")
  }
  
  # Make sure the selected scale and subcompartments exist
  if(!scale %in% unique(solution$Scale)){
    stop("Selected scale does not exist")
  }
  
  # Aggregate over species
  cnames <- names(solution)
  cnames <- cnames[!cnames %in% c("Species", "Mass_kg")]
  formula <- as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + ")))
  solution <- aggregate(formula, data = solution, sum)
  
  if (!is.null(scale)) {
    solution <- solution[solution$Scale %in% scale, ]
  } else {
    solution <- solution
  }
  
  solution <- solution |>
    group_by(SubCompart, Scale) |>
    summarise(Mass_kg = mean(Mass_kg)) |>
    ungroup()
  
  ggplot(solution, aes(area = Mass_kg, fill = Mass_kg, label = SubCompart)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white",
                                  place = "centre",
                                  size = 15) +
    labs(title = paste0("Distribution of steady state mean masses at ", scale, " scale"),
         fill = "Mass [kg]") 
}