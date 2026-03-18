#####################
# Plotting functions
#####################

# Lookup table for legend labels
subcompart_labels <- c(
  agriculturalsoil    = "Agricultural soil",
  air                 = "Air",
  cloudwater          = "Cloud water",
  freshwatersediment  = "Freshwater sediment",
  lake                = "Lake",
  marinesediment      = "Marine sediment",
  naturalsoil         = "Natural soil",
  othersoil           = "Other soil",
  river               = "River",
  sea                 = "Sea",
  deepocean           = "Deep ocean",
  lakesediment        = "Lake sediment"
)

# Subcompartment colors
subcompart_colors <- c(
  # Marine waters
  "Sea"        = "#FF69B4",  
  "Deep ocean" = "#8A2BE2",  
  
  # Freshwaters
  "River"      = "#87CEEB",
  "Lake"       = "#4A90E2",  
  
  # Atmosphere
  "Air"                 = "#b2df8a",
  "Cloud water"         = "#006400",
  
  # Soils
  "Agricultural soil"   = "#d9a066",
  "Natural soil"        = "#8c510a",
  "Other soil"          = "#3b2a1a",
  
  # Sediments
  "Freshwater sediment" = "#ffe066",
  "Marine sediment"     = "#fd8d3c",
  "Lake sediment"       = "#CC0000"
)


################################################################################################################################################
################################################################################################################################################
# Solution plots
################################################################################################################################################
################################################################################################################################################

####################################################################
# Deterministic & steady state
####################################################################
DetSSPlot <- function(scale = NULL, subcompart = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) 
    stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(solution), c("Species","Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse=" + "))), data=solution, sum)
  
  solution <- solution[solution$Scale == scale, ]
  if(!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution, aes(x=SubCompartLabel, y=Mass_kg, fill=SubCompartLabel)) +
    geom_col() +
    theme_bw() +
    labs(title=paste0("Steady state mass at ", scale, " scale"),
         x="", y=paste0("Mass of ", World$substance, " [kg]")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n=10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_manual(values=subcompart_colors) +
    theme(legend.position="right",
          axis.text.x = element_text(angle=45, hjust=1))
}

####################################################################
# Deterministic & dynamic
####################################################################
DetDynSolPlot <- function(scale = NULL, subcompart = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'time', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(solution), c("Species","Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse=" + "))), data=solution, sum)
  
  solution <- solution[solution$Scale==scale, ]
  if(!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  solution$Year <- solution$time/(365.25*24*3600)
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution, aes(x=Year, y=Mass_kg, group=SubCompartLabel, color=SubCompartLabel)) +
    geom_line() +
    theme_bw() +
    labs(title=paste0("Dynamic mass at ", scale, " scale"),
         x="Year", y=paste0("Mass of ", World$substance, " [kg]")) +
    scale_color_manual(values=subcompart_colors) +
    guides(color=guide_legend(title="Subcompartment")) +
    scale_y_continuous(labels=scales::label_scientific())
}

####################################################################
# Probabilistic & dynamic
####################################################################
ProbDynSolPlot <- function(scale=NULL, subcompart=NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  solution <- solution[c('SubCompart','Scale','Species','time','RUNs','Mass_kg')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) 
    stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(solution), c("Species","Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse=" + "))), data=solution, sum)
  
  solution <- solution[solution$Scale==scale, ]
  if(!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  solution$time <- as.numeric(as.character(solution$time))
  solution$Year <- solution$time/(365.25*24*3600)
  
  summary_stats <- solution |>
    group_by(Year, SubCompart) |>
    summarise(
      Mean_Value = mean(Mass_kg, na.rm=TRUE),
      SD_Value   = sd(Mass_kg, na.rm=TRUE)
    ) |>
    ungroup() |>
    mutate(
      Lower_CI = Mean_Value - 1.96*SD_Value/sqrt(n()),
      Upper_CI = Mean_Value + 1.96*SD_Value/sqrt(n()),
      SubCompartLabel = subcompart_labels[SubCompart]
    )
  
  ggplot(summary_stats, aes(x=Year, y=Mean_Value, color=SubCompartLabel, fill=SubCompartLabel)) +
    geom_ribbon(aes(ymin=Lower_CI, ymax=Upper_CI), alpha=0.2, show.legend=FALSE) +
    geom_line(size=1) +
    theme_minimal() +
    labs(title=paste0("Dynamic mean mass at ", scale, " scale"),
         subtitle="with uncertainty bands over time",
         x="Year",
         y=paste0("Mass of ", World$substance, " [kg]")) +
    scale_color_manual(values=subcompart_colors) +
    scale_fill_manual(values=subcompart_colors) +
    guides(color=guide_legend(title="Subcompartment"))
}

####################################################################
# Probabilistic & steady state
####################################################################
ProbSSSolPlot <- function(scale=NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  solution <- solution[c('SubCompart','Scale','Species','RUNs','Mass_kg')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  nRUNs <- length(unique(solution$RUNs))
  
  solution <- solution |>
    group_by(SubCompart, Scale, RUNs) |>
    summarise(Mass_kg = sum(Mass_kg)) |>
    ungroup() |>
    group_by(SubCompart, Scale) |>
    summarise(Mass_kg = mean(Mass_kg), n=n()) |>
    ungroup()
  
  if(nRUNs != unique(solution$n)) stop("nRUNs not equal to n in summarise")
  
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution[solution$Scale==scale,], aes(area=Mass_kg, fill=SubCompartLabel,
                                               label=paste(SubCompartLabel, round(Mass_kg/sum(Mass_kg)*100,2), "%", sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_manual(values=subcompart_colors) +
    labs(title=paste0("Distribution of average steady state masses at ", scale, " scale")) +
    theme(legend.position="right")
}


################################################################################################################################################
################################################################################################################################################
# Concentration plots
################################################################################################################################################
################################################################################################################################################

####################################################################
# Deterministic & steady state
####################################################################
DetSSConcPlot <- function(scale=NULL, subcompart=NULL){
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','Concentration','Unit')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(conc$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  conc <- conc[conc$Scale==scale, ]
  if(!is.null(subcompart)) conc <- conc[conc$SubCompart %in% subcompart, ]
  
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  ggplot(conc, aes(x=SubCompartLabel, y=Concentration, fill=SubCompartLabel)) +
    geom_col() +
    theme_bw() +
    labs(title=paste0("Steady state concentration at ", scale, " scale"),
         x="", y=paste0("Concentration of ", World$substance)) +
    scale_fill_manual(values=subcompart_colors) +
    theme(legend.position="right",
          axis.text.x=element_text(angle=45, hjust=1))
}

####################################################################
# Deterministic & dynamic
####################################################################
DetDynConcPlot <- function(scale=NULL, subcompart=NULL){
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','time','Concentration','Unit')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(conc$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  conc <- conc[conc$Scale==scale, ]
  if(!is.null(subcompart)) conc <- conc[conc$SubCompart %in% subcompart, ]
  
  conc$Year <- conc$time/(365.25*24*3600)
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  ggplot(conc, aes(x=Year, y=Concentration, group=SubCompartLabel, color=SubCompartLabel)) +
    geom_line() +
    theme_bw() +
    labs(title=paste0("Dynamic concentration at ", scale, " scale"),
         x="Year", y=paste0("Concentration of ", World$substance)) +
    scale_color_manual(values=subcompart_colors) +
    guides(color=guide_legend(title="Subcompartment"))
}

####################################################################
# Probabilistic & dynamic
####################################################################
ProbDynConcPlot <- function(scale = NULL, subcompart = NULL){
  
  # Merge concentration with state info
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','time','RUNs','Concentration','Unit')]
  
  # Validation
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(conc$SubCompart))) 
    stop("One or more selected subcompartments do not exist")
  
  # Aggregate over species
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  # Filter for chosen scale and subcompartments
  conc <- conc[conc$Scale == scale, ]
  if(!is.null(subcompart)) conc <- conc[conc$SubCompart %in% subcompart, ]
  
  # Convert time to numeric & to years
  conc$time <- as.numeric(as.character(conc$time))
  conc$Year <- conc$time / (365.25*24*3600)
  
  # Compute mean and SD for each subcompartment per year
  summary_stats <- conc |>
    group_by(Year, SubCompart) |>
    summarise(
      Mean_Value = mean(Concentration, na.rm=TRUE),
      SD_Value   = sd(Concentration, na.rm=TRUE),
      .groups="drop"
    ) |>
    mutate(
      Lower_CI = Mean_Value - 1.96*SD_Value/sqrt(n()),
      Upper_CI = Mean_Value + 1.96*SD_Value/sqrt(n()),
      SubCompartLabel = subcompart_labels[SubCompart]
    )
  
  ggplot(summary_stats, aes(x=Year, y=Mean_Value, color=SubCompartLabel, fill=SubCompartLabel)) +
    geom_ribbon(aes(ymin=Lower_CI, ymax=Upper_CI), alpha=0.2, show.legend=FALSE) +
    geom_line(size=1) +
    theme_minimal() +
    labs(
      title = paste0("Dynamic mean concentration at ", scale, " scale"),
      subtitle = "with uncertainty bands over time",
      x = "Year",
      y = paste0("Concentration of ", World$substance, " [", unique(conc$Unit), "]")
    ) +
    scale_color_manual(values = subcompart_colors) +
    scale_fill_manual(values = subcompart_colors) +
    guides(color = guide_legend(title="Subcompartment"))
}

####################################################################
# Probabilistic & steady state
####################################################################
ProbSSConcPlot <- function(scale=NULL){
  
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','RUNs','Concentration','Unit')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  
  # Aggregate over species and runs
  nRUNs <- length(unique(conc$RUNs))
  
  conc <- conc |>
    group_by(SubCompart, Scale, RUNs) |>
    summarise(Concentration = sum(Concentration), .groups="drop") |>
    group_by(SubCompart, Scale) |>
    summarise(Concentration = mean(Concentration), n = n(), .groups="drop")
  
  if(nRUNs != unique(conc$n)) stop("nRUNs not equal to n in summarise")
  
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  ggplot(conc[conc$Scale==scale,], aes(x=SubCompartLabel, y=Concentration, fill=SubCompartLabel,
                                       label=paste0(SubCompartLabel, " (", round(Concentration/sum(Concentration)*100,2),"%)"))) +
    geom_violin() +
    theme_bw() +
    labs(title=paste0("Average steady state concentration at ", scale, " scale"),
         x="", y=paste0("Concentration of ", World$substance)) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n=10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_manual(values=subcompart_colors) +
    theme(legend.position="right",
          axis.text.x=element_text(angle=45, hjust=1))
}


################################################################################################################################################
################################################################################################################################################
# Mass distribution plots
################################################################################################################################################
################################################################################################################################################

####################################################################
# Deterministic & steady state
####################################################################
DetSSMassDist <- function(scale = NULL){
  
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  solution <- solution[c("SubCompart","Scale","Species","Mass_kg")]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  # Aggregate over species
  cnames <- setdiff(names(solution), c("Species","Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse=" + "))), data=solution, sum)
  
  # Treemap over scales
  scale_solution <- solution |>
    group_by(Scale) |>
    summarise(Mass_kg = sum(Mass_kg), .groups="drop") |>
    mutate(Mass_percent = round((Mass_kg/sum(Mass_kg))*100,2),
           Mass_percent_label = paste0(Mass_percent,"%"))
  
  scale_plot <- ggplot(scale_solution, aes(area=Mass_kg, fill=Scale,
                                           label=paste(Scale, Mass_percent_label, sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_discrete() +
    labs(title="Distribution of steady state masses over scales") +
    theme(legend.position="none")
  
  # Treemap for subcompartments at chosen scale
  subcompart_solution <- solution[solution$Scale==scale, ] |>
    mutate(Mass_percent = round((Mass_kg/sum(Mass_kg))*100,2),
           Mass_percent_label = paste0(Mass_percent,"%"),
           SubCompartLabel = subcompart_labels[SubCompart])
  
  subcompart_plot <- ggplot(subcompart_solution, aes(area=Mass_kg, fill=SubCompartLabel,
                                                     label=paste(SubCompartLabel, Mass_percent_label, sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_manual(values=subcompart_colors) +
    labs(title=paste0("Distribution of steady state masses at ", scale, " scale")) +
    theme(legend.position="right")
  
  grid.arrange(scale_plot, subcompart_plot, ncol=1)
}

####################################################################
# Probabilistic & steady state
####################################################################
ProbSSMassDist <- function(scale = NULL){
  
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  nRUNs <- length(unique(solution$RUNs))
  
  solution <- solution |>
    group_by(SubCompart, Scale, RUNs) |>
    summarise(Mass_kg = sum(Mass_kg), .groups="drop") |>
    group_by(SubCompart, Scale) |>
    summarise(Mass_kg = mean(Mass_kg), n=n(), .groups="drop")
  
  if(nRUNs != unique(solution$n)) stop("nRUNs not equal to n in summarise")
  
  # Top-level treemap over scales
  scale_solution <- solution |>
    group_by(Scale) |>
    summarise(Mass_kg = sum(Mass_kg), .groups="drop") |>
    mutate(Mass_percent = round((Mass_kg/sum(Mass_kg))*100,2),
           Mass_percent_label = paste0(Mass_percent,"%"))
  
  scale_plot <- ggplot(scale_solution, aes(area=Mass_kg, fill=Scale,
                                           label=paste(Scale, Mass_percent_label, sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_discrete() +
    labs(title="Distribution of average steady state masses over scales") +
    theme(legend.position="none")
  
  # Subcompartment treemap
  subcompart_solution <- solution[solution$Scale==scale, ] |>
    mutate(Mass_percent = round((Mass_kg/sum(Mass_kg))*100,2),
           Mass_percent_label = paste0(Mass_percent,"%"),
           SubCompartLabel = subcompart_labels[SubCompart])
  
  subcompart_plot <- ggplot(subcompart_solution, aes(area=Mass_kg, fill=SubCompartLabel,
                                                     label=paste(SubCompartLabel, Mass_percent_label, sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_manual(values=subcompart_colors) +
    labs(title=paste0("Distribution of average steady state masses at ", scale, " scale")) +
    theme(legend.position="right")
  
  grid.arrange(scale_plot, subcompart_plot, ncol=1)
}






