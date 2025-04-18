% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6SBcore.R
\name{SBcore}
\alias{SBcore}
\title{SBcore}
\description{
the centre object for running simple box with an R script.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{states}}{getter for r.o. property}

\item{\code{kaas}}{getter for r.o. property (all k's)}

\item{\code{nodelist}}{getter for r.o. property}

\item{\code{moduleList}}{getter for r.o. property}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SBcore-new}{\code{SBcore$new()}}
\item \href{#method-SBcore-NewProcess}{\code{SBcore$NewProcess()}}
\item \href{#method-SBcore-PostponeVarProcess}{\code{SBcore$PostponeVarProcess()}}
\item \href{#method-SBcore-NewFlow}{\code{SBcore$NewFlow()}}
\item \href{#method-SBcore-NewCalcVariable}{\code{SBcore$NewCalcVariable()}}
\item \href{#method-SBcore-NewSolver}{\code{SBcore$NewSolver()}}
\item \href{#method-SBcore-FindStatefrom3D}{\code{SBcore$FindStatefrom3D()}}
\item \href{#method-SBcore-Solve}{\code{SBcore$Solve()}}
\item \href{#method-SBcore-makeInvFuns}{\code{SBcore$makeInvFuns()}}
\item \href{#method-SBcore-exportEngine}{\code{SBcore$exportEngine()}}
\item \href{#method-SBcore-exportEngineR}{\code{SBcore$exportEngineR()}}
\item \href{#method-SBcore-Masses}{\code{SBcore$Masses()}}
\item \href{#method-SBcore-VariableValues}{\code{SBcore$VariableValues()}}
\item \href{#method-SBcore-Concentration}{\code{SBcore$Concentration()}}
\item \href{#method-SBcore-Emissions}{\code{SBcore$Emissions()}}
\item \href{#method-SBcore-PlotConcentration}{\code{SBcore$PlotConcentration()}}
\item \href{#method-SBcore-PlotMasses}{\code{SBcore$PlotMasses()}}
\item \href{#method-SBcore-PlotMassDistribution}{\code{SBcore$PlotMassDistribution()}}
\item \href{#method-SBcore-MassesAsRelational}{\code{SBcore$MassesAsRelational()}}
\item \href{#method-SBcore-metaData}{\code{SBcore$metaData()}}
\item \href{#method-SBcore-fetchDims}{\code{SBcore$fetchDims()}}
\item \href{#method-SBcore-fetchData}{\code{SBcore$fetchData()}}
\item \href{#method-SBcore-fetch_current}{\code{SBcore$fetch_current()}}
\item \href{#method-SBcore-fetchDataUnits}{\code{SBcore$fetchDataUnits()}}
\item \href{#method-SBcore-doInherit}{\code{SBcore$doInherit()}}
\item \href{#method-SBcore-findState}{\code{SBcore$findState()}}
\item \href{#method-SBcore-filterStatesFrame}{\code{SBcore$filterStatesFrame()}}
\item \href{#method-SBcore-whichDataTable}{\code{SBcore$whichDataTable()}}
\item \href{#method-SBcore-UpdateKaas}{\code{SBcore$UpdateKaas()}}
\item \href{#method-SBcore-VarsFromprocesses}{\code{SBcore$VarsFromprocesses()}}
\item \href{#method-SBcore-CalcVar}{\code{SBcore$CalcVar()}}
\item \href{#method-SBcore-SetConst}{\code{SBcore$SetConst()}}
\item \href{#method-SBcore-UpdateDirty}{\code{SBcore$UpdateDirty()}}
\item \href{#method-SBcore-DependOn}{\code{SBcore$DependOn()}}
\item \href{#method-SBcore-whichUnresolved}{\code{SBcore$whichUnresolved()}}
\item \href{#method-SBcore-CleanupCalcGraphAbove}{\code{SBcore$CleanupCalcGraphAbove()}}
\item \href{#method-SBcore-mutateVars}{\code{SBcore$mutateVars()}}
\item \href{#method-SBcore-UpdateData}{\code{SBcore$UpdateData()}}
\item \href{#method-SBcore-FromDataAndTo}{\code{SBcore$FromDataAndTo()}}
\item \href{#method-SBcore-save_world}{\code{SBcore$save_world()}}
\item \href{#method-SBcore-load_world}{\code{SBcore$load_world()}}
\item \href{#method-SBcore-clone}{\code{SBcore$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-new"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-new}{}}}
\subsection{Method \code{new()}}{
init
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$new(NewstateModule)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{NewstateModule}}{The stateModule with its States and standard input data}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-NewProcess"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-NewProcess}{}}}
\subsection{Method \code{NewProcess()}}{
add a process to the calculations
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$NewProcess(ProcessFunction)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ProcessFunction}}{The name (as character) of the process defining function 
see}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-PostponeVarProcess"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-PostponeVarProcess}{}}}
\subsection{Method \code{PostponeVarProcess()}}{
remove a set of SB items from the calculations, for later reprocessing
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$PostponeVarProcess(
  VarFunctions = NULL,
  FlowFunctions = NULL,
  ProcesFunctions
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{VarFunctions}}{The name (as character) of the process defining function 2be postponed}

\item{\code{FlowFunctions}}{The name (as character) of the process defining function 2be postponed}

\item{\code{ProcesFunctions}}{The name (as character) of the process defining function 2be postponed}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The SBitems will be removed from the internal structures, but the names of their defining functions 
will be stored locally. The private method DoPostponed will add and execute them. (This is automised, 
because it will be needed in sensitivity analyses etc.)
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-NewFlow"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-NewFlow}{}}}
\subsection{Method \code{NewFlow()}}{
add a flow to the calculations
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$NewFlow(FlowFunction, WithProcess = "k_Advection")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{FlowFunction}}{The name (as character) of the flow defining function 
see}

\item{\code{WithProcess}}{The name (as character) of the process related to this function
This is normally the k_Advection process; created if not already present}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-NewCalcVariable"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-NewCalcVariable}{}}}
\subsection{Method \code{NewCalcVariable()}}{
add a SBVariable to the calculations
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$NewCalcVariable(VariableFunction, AggrBy = NA, AggrFun = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{VariableFunction}}{The name (as character) of the SBvariable defining function 
see}

\item{\code{AggrBy}}{exceptional parameter to aggregate over a dim (1 of The3D)}

\item{\code{AggrFun}}{FUN to use in aggregation over a dim (1 of The3D)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-NewSolver"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-NewSolver}{}}}
\subsection{Method \code{NewSolver()}}{
create and set the solver (there can only be one)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$NewSolver(SolverFunction, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{SolverFunction}}{The (name of) the solver defining function}

\item{\code{...}}{passed on to init of SolverModule i.e. to CalcGraphModule$MoreParams
for use in the defining function}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-FindStatefrom3D"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-FindStatefrom3D}{}}}
\subsection{Method \code{FindStatefrom3D()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$FindStatefrom3D(df3Ds)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-Solve"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-Solve}{}}}
\subsection{Method \code{Solve()}}{
Run the matrix with exchange constants; 
the exact calculation is described in the function that defined the solve module.
This can be solving for the steady state or running the system for a period of time.
Results are available
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$Solve(emissions, needdebug = F, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{emissions}}{data.frame-ish with columns "Abbr" (state defining character) and 
"emis" numbers}

\item{\code{needdebug}}{if T the defining function will open in debugging modus}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-makeInvFuns"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-makeInvFuns}{}}}
\subsection{Method \code{makeInvFuns()}}{
Creates inverse functions for each row of the given dataframe
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$makeInvFuns(paramdf)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{paramdf}}{A dataframe with one row for eacht distribution. Must at 
least contain the columns "Distribution" (which can contain "triangular", "uniform" or "normal"), 
"a", "b", and "c".

"a", "b" and "c" contain the following:
triangular: a = Minimum, b = Maximum, c = Peak value
normal: a = Mean, b = Standard deviation 
uniform: a = Minimum , b = Maximum}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-exportEngine"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-exportEngine}{}}}
\subsection{Method \code{exportEngine()}}{
Export the matrix of speed constants, aka Engine, to an excel file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$exportEngine(excelFile)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-exportEngineR"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-exportEngineR}{}}}
\subsection{Method \code{exportEngineR()}}{
Export the matrix to World
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$exportEngineR()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-Masses"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-Masses}{}}}
\subsection{Method \code{Masses()}}{
Save the last calculated masses in the core
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$Masses()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-VariableValues"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-VariableValues}{}}}
\subsection{Method \code{VariableValues()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$VariableValues()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-Concentration"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-Concentration}{}}}
\subsection{Method \code{Concentration()}}{
Function to obtain steady state concentrations, using the solution saved in world.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$Concentration()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-Emissions"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-Emissions}{}}}
\subsection{Method \code{Emissions()}}{
Save the last used emissions in the core
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$Emissions()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-PlotConcentration"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-PlotConcentration}{}}}
\subsection{Method \code{PlotConcentration()}}{
Return the appropriate concentration plot
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$PlotConcentration(scale = NULL, subcompart = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-PlotMasses"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-PlotMasses}{}}}
\subsection{Method \code{PlotMasses()}}{
Return the appropriate solution plot
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$PlotMasses(scale = NULL, subcompart = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-PlotMassDistribution"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-PlotMassDistribution}{}}}
\subsection{Method \code{PlotMassDistribution()}}{
Return the appropriate mass distribution plot
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$PlotMassDistribution(scale = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-MassesAsRelational"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-MassesAsRelational}{}}}
\subsection{Method \code{MassesAsRelational()}}{
Injection from SolverModule
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$MassesAsRelational(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-metaData"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-metaData}{}}}
\subsection{Method \code{metaData()}}{
Obtain the names of the variables and tablename in which the data resides 
Not needed in normal use
see
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$metaData()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-fetchDims"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-fetchDims}{}}}
\subsection{Method \code{fetchDims()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$fetchDims(vars)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-fetchData"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-fetchData}{}}}
\subsection{Method \code{fetchData()}}{
Obtain the data for a SBvariable or a flow
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$fetchData(varname = "all")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{varname}}{the name of the variable. Returns a list of variables if varname == "all"}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-fetch_current"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-fetch_current}{}}}
\subsection{Method \code{fetch_current()}}{
fetch specific values from core
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$fetch_current(withoutValues)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{withoutValues}}{data.frame-ish with columns `varname` and needed D}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-fetchDataUnits"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-fetchDataUnits}{}}}
\subsection{Method \code{fetchDataUnits()}}{
function to obtain the data for a variable or flow, including the units whenever present in the Units csv
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$fetchDataUnits(varname = "all")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{varname}}{name of the variable}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-doInherit"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-doInherit}{}}}
\subsection{Method \code{doInherit()}}{
Pseudo-inherit data; for instance from compartment to subcompartments.
Not needed in normal use.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$doInherit(fromData, toData)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fromData}}{the name of the variable before copying}

\item{\code{toData}}{the name of the variable after copying}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-findState"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-findState}{}}}
\subsection{Method \code{findState()}}{
returns the state for a state-name. Injected from States property
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$findState(abbr)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{abbr}}{the name of the state}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-filterStatesFrame"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-filterStatesFrame}{}}}
\subsection{Method \code{filterStatesFrame()}}{
limits the states according to the filter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$filterStatesFrame(inDataFrame)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{the name of one of the dimensions of the states (ScaleName, SubCompartName, SpeciesName)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-whichDataTable"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-whichDataTable}{}}}
\subsection{Method \code{whichDataTable()}}{
returns the table, determined by the dimensions et al.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$whichDataTable(KeyNames)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{KeyNames}}{the dim (of The3D etc) determine which table}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-UpdateKaas"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-UpdateKaas}{}}}
\subsection{Method \code{UpdateKaas()}}{
runs (or tries to) the calculation graph, either the whole graph or 
limited to a single process
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$UpdateKaas(aProcessModule = NULL, mergeExisting = T)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{aProcessModule}}{the name of the state}

\item{\code{mergeExisting}}{normally leaves all other processes unchanged, but if 
not mergeExisting, all process results (kaas) are cleared first}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-VarsFromprocesses"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-VarsFromprocesses}{}}}
\subsection{Method \code{VarsFromprocesses()}}{
Tries to create all variable modules that current processes and all flows need
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$VarsFromprocesses()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-CalcVar"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-CalcVar}{}}}
\subsection{Method \code{CalcVar()}}{
runs (or tries to) the calculation for a Variable and stores the results.
After this, the SBvariable can be viewed with fetchData and this data will be used 
when calculating processes, or other variables.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$CalcVar(aVariable)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{aVariable}}{the name of the Variable}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-SetConst"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-SetConst}{}}}
\subsection{Method \code{SetConst()}}{
set a constant in the internal data, to enable use by SB variable etc.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$SetConst(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{named value}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-UpdateDirty"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-UpdateDirty}{}}}
\subsection{Method \code{UpdateDirty()}}{
runs (or tries to) the calculation for Variables,
and continues from there to update all processes and variables 
that have a depency of any of the Variables, recursively.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$UpdateDirty(Variables)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{Variables}}{the name(s) of the Variable(s) that are "dirty" 
(Datalayer has been updated)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-DependOn"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-DependOn}{}}}
\subsection{Method \code{DependOn()}}{
Which Graph elements SBVars depend on?
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$DependOn(SBVars)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-whichUnresolved"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-whichUnresolved}{}}}
\subsection{Method \code{whichUnresolved()}}{
Verifies the presence of needed variables for the calculation of 
all processes
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$whichUnresolved()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-CleanupCalcGraphAbove"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-CleanupCalcGraphAbove}{}}}
\subsection{Method \code{CleanupCalcGraphAbove()}}{
removes variables appearently not needed in the calculation DAG
starting from Varname upwards
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$CleanupCalcGraphAbove(VarName)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{VarName}}{the name of the Variable which inputs are cleaned, recursively}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-mutateVars"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-mutateVars}{}}}
\subsection{Method \code{mutateVars()}}{
if UpdateRows is a list of variables (each containing a fetchdata() result): see private$MutateVar,
if UpdateRows is (a csv filename of) a single dataframe it will be converted before the call
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$mutateVars(UpdateRows)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-UpdateData"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-UpdateData}{}}}
\subsection{Method \code{UpdateData()}}{
Replaces a complete table in the internal data system. 
Use with care.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$UpdateData(UpdateDF, keys, TableName = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{UpdateDF}}{the new table, replcing the previous}

\item{\code{keys}}{dims that determines the table to be updated, see whichDataTable}

\item{\code{TableName}}{alternative to keys: provide the tablename. Use with care!}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-FromDataAndTo"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-FromDataAndTo}{}}}
\subsection{Method \code{FromDataAndTo()}}{
derive transfer processes from datalayer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$FromDataAndTo(processName = "all")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{processName}}{the one you are looking for, or (default) for "all" processes}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-save_world"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-save_world}{}}}
\subsection{Method \code{save_world()}}{
states and all data-layer to an rds-file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$save_world(filename)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-load_world"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-load_world}{}}}
\subsection{Method \code{load_world()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$load_world(filename)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SBcore-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SBcore-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SBcore$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
