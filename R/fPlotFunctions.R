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


subcompart_colors <- c(
  # Waters
  "Sea"                 = "#1f78b4",   # blue
  "River"               = "#a6cee3",   # light blue
  "Lake"                = "#3690c0",   # medium blue
  "Deep ocean"          = "#08306b",   # dark blue
  
  # Atmosphere
  "Air"                 = "#b2df8a",   # light green
  "Cloud water"         = "#33a02c",   # medium green
  
  # Soils
  "Agricultural soil"   = "#c18750",   # light brown
  "Natural soil"        = "#8b5a2b",   # medium brown
  "Other soil"          = "#5c3317",   # dark brown
  
  # Sediments
  "Freshwater sediment" = "#ffffb2",   # pale yellow
  "Marine sediment"     = "#fecc5c",   # medium yellow
  "Lake sediment"       = "#fd8d3c"    # darker yellow/orange
)



############################################################################################################################################################
############################################################################################################################################################
# Functions for solution plots 
############################################################################################################################################################
############################################################################################################################################################

# Deterministic & steady state
DetSSPlot <- function(scale = NULL, subcompart = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  # Aggregate over species
  cnames <- setdiff(names(solution), c("Species", "Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + "))), data = solution, sum)
  
  solution <- solution[solution$Scale == scale, ]
  if (!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  # Map labels
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution, aes(x = SubCompartLabel, y = Mass_kg, fill = SubCompartLabel)) +
    geom_col() +
    theme_bw() +
    labs(title = paste0("Steady state mass at ", scale, " scale"),
         x = "", y = paste0("Mass of ", World$substance, " [kg]")) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_manual(values = subcompart_colors) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deterministic & dynamic
DetDynSolPlot <- function(scale = NULL, subcompart = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'time', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(solution), c("Species", "Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + "))), data = solution, sum)
  
  solution <- solution[solution$Scale == scale, ]
  if (!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  # Convert time
  solution$Year <- solution$time / (365.25*24*3600)
  
  # Map labels
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution, aes(x = Year, y = Mass_kg, group = SubCompartLabel, color = SubCompartLabel)) +
    geom_line() +
    theme_bw() +
    labs(title = paste0("Dynamic mass at ", scale, " scale"),
         x = "Year", y = paste0("Mass of ", World$substance, " [kg]")) +
    scale_y_continuous(labels = scales::label_scientific()) +
    scale_color_manual(values = subcompart_colors) +
    guides(color = guide_legend(title = "Subcompartment"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Probabilistic & dynamic
ProbDynSolPlot <- function(scale = NULL, subcompart = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'time', 'RUNs', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  if (!is.null(subcompart) && !all(subcompart %in% unique(solution$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(solution), c("Species", "Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse = " + "))), data = solution, sum)
  
  solution <- solution[solution$Scale == scale, ]
  if (!is.null(subcompart)) solution <- solution[solution$SubCompart %in% subcompart, ]
  
  # Convert time
  solution$Year <- solution$time / (365.25*24*3600)
  
  # Summary stats
  summary_stats <- solution |>
    group_by(Year, SubCompart) |>
    summarise(Mean_Value = mean(Mass_kg), SD_Value = sd(Mass_kg)) |>
    ungroup() |>
    mutate(Lower_CI = Mean_Value - 1.96*SD_Value/sqrt(n()),
           Upper_CI = Mean_Value + 1.96*SD_Value/sqrt(n()))
  
  summary_stats$SubCompartLabel <- subcompart_labels[summary_stats$SubCompart]
  
  ggplot(summary_stats, aes(x = Year, y = Mean_Value, color = SubCompartLabel)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = SubCompartLabel), alpha = 0.2, show.legend = FALSE) +
    theme_minimal() +
    labs(title = paste0("Dynamic mean mass in ", paste(subcompart, collapse=", "), " at ", scale, " scale"),
         subtitle = "with uncertainty bands over time",
         x = "Year", y = paste0("Mass of ", World$substance, " [kg]")) +
    scale_color_manual(values = subcompart_colors) +
    guides(color = guide_legend(title = "Subcompartment"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Probabilistic & steady state
ProbSSSolPlot <- function(scale = NULL){
  solution <- merge(World$Masses(), World$states$asDataFrame, by = "Abbr")
  solution <- solution[c('SubCompart', 'Scale', 'Species', 'RUNs', 'Mass_kg')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  nRUNs <- length(unique(solution$RUNs))
  
  solution <- solution |>
    group_by(SubCompart, Scale, RUNs) |>
    summarise(Mass_kg = sum(Mass_kg)) |>
    ungroup() |>
    group_by(SubCompart, Scale) |>
    summarise(Mass_kg = mean(Mass_kg), n = n()) |>
    ungroup()
  
  if(nRUNs != unique(solution$n)) stop("nRUNs not equal to n in summarise")
  
  solution$SubCompartLabel <- subcompart_labels[solution$SubCompart]
  
  ggplot(solution[solution$Scale==scale,], aes(area = Mass_kg, fill = SubCompartLabel,
                                               label = paste(SubCompartLabel, round(Mass_kg/sum(Mass_kg)*100,2), "%", sep="\n"))) +
    geom_treemap() +
    geom_treemap_text(colour="white", place="centre", size=15) +
    scale_fill_manual(values=subcompart_colors) +
    labs(title = paste0("Distribution of average steady state masses at ", scale, " scale")) +
    theme(legend.position="right")
}



############################################################################################################################################################
############################################################################################################################################################
# Functions for concentration plots
############################################################################################################################################################
############################################################################################################################################################

# Deterministic & steady state
DetSSConcPlot <- function(scale = NULL, subcompart = NULL){
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','Concentration','Unit')]
  
  if(length(scale) != 1) stop("Please select 1 scale")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(conc$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  conc <- conc[conc$Scale == scale, ]
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deterministic & dynamic
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
  
  conc$Year <- conc$time / (365.25*24*3600)
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  ggplot(conc, aes(x=Year, y=Concentration, group=SubCompartLabel, color=SubCompartLabel)) +
    geom_line() +
    theme_bw() +
    labs(title=paste0("Dynamic concentration at ", scale, " scale"),
         x="Year", y=paste0("Concentration of ", World$substance)) +
    scale_y_continuous(labels=scales::label_scientific()) +
    scale_color_manual(values=subcompart_colors) +
    guides(color=guide_legend(title="Subcompartment"))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Probabilistic & dynamic
ProbDynConcPlot <- function(scale=NULL, subcompart=NULL){
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','time','RUNs','Concentration','Unit')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(length(subcompart)!=1) stop("Please select 1 subcompartment")
  if(!scale %in% unique(conc$Scale)) stop("Selected scale does not exist")
  if(!is.null(subcompart) && !all(subcompart %in% unique(conc$SubCompart))) stop("One or more selected subcomparts do not exist")
  
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  conc <- conc[conc$Scale==scale, ]
  if(!is.null(subcompart)) conc <- conc[conc$SubCompart %in% subcompart, ]
  
  conc$Year <- conc$time / (365.25*24*3600)
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  summary_stats <- conc |>
    group_by(Year, SubCompartLabel) |>
    summarise(
      Mean_Value = mean(Concentration, na.rm=TRUE),
      SD_Value = sd(Concentration, na.rm=TRUE),
      .groups="drop"
    ) |>
    mutate(Lower_CI = Mean_Value - 1.96*SD_Value/sqrt(n()),
           Upper_CI = Mean_Value + 1.96*SD_Value/sqrt(n()))
  
  ggplot(summary_stats, aes(x=Year, y=Mean_Value, color=SubCompartLabel)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=Lower_CI, ymax=Upper_CI, fill=SubCompartLabel), alpha=0.2, show.legend=FALSE) +
    theme_minimal() +
    labs(title=paste0("Dynamic mean concentration in ", subcompart, " at ", scale, " scale"),
         subtitle="with uncertainty bands over time",
         x="Year", y=paste0("Concentration of ", World$substance)) +
    scale_color_manual(values=subcompart_colors) +
    guides(color=guide_legend(title="Subcompartment"))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Probabilistic & steady state 
ProbSSConcPlot <- function(scale=NULL){
  conc <- merge(World$Concentration(), World$states$asDataFrame, by="Abbr")
  conc <- conc[c('SubCompart','Scale','Species','RUNs','Concentration','Unit')]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  
  cnames <- setdiff(names(conc), c("Species","Concentration"))
  conc <- aggregate(as.formula(paste("Concentration ~", paste(cnames, collapse=" + "))), data=conc, sum)
  
  conc$SubCompartLabel <- subcompart_labels[conc$SubCompart]
  
  ggplot(conc[conc$Scale==scale,], aes(x=SubCompartLabel, y=Concentration, fill=SubCompartLabel,
                                       label=paste0(SubCompartLabel, " (", round(Concentration/sum(Concentration)*100,2),"%)"))) +
    geom_violin() +
    theme_bw() +
    labs(title=paste0("Steady state concentration at ", scale, " scale"),
         x="", y=paste0("Concentration of ", World$substance)) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 10),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_manual(values=subcompart_colors) +
    theme(legend.position="right",
          axis.text.x=element_text(angle=45, hjust=1))
}



############################################################################################################################################################
############################################################################################################################################################
# Functions for mass distribution plots
############################################################################################################################################################
############################################################################################################################################################

# Deterministic & steady state
DetSSMassDist <- function(scale = NULL){
  
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  solution <- solution[c("SubCompart","Scale","Species","Mass_kg")]
  
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  # Aggregate over species
  cnames <- setdiff(names(solution), c("Species","Mass_kg"))
  solution <- aggregate(as.formula(paste("Mass_kg ~", paste(cnames, collapse=" + "))), data=solution, sum)
  
  # Aggregate over scales for top-level treemap
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
  
  # Filter for selected scale
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


# Probabilistic & steady state
ProbSSMassDist <- function(scale = NULL){
  
  solution <- merge(World$Masses(), World$states$asDataFrame, by="Abbr")
  if(length(scale)!=1) stop("Please select 1 scale")
  if(!scale %in% unique(solution$Scale)) stop("Selected scale does not exist")
  
  # Aggregate over species and average over runs
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
  
  # Filter for selected scale and prepare subcompartment treemap
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

