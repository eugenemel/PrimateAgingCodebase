#pipeline
LIBRARIES=c("readxl","dplyr","survival","ggplot2","survminer","thematic","flexsurv","ggrepel","phytools","brms","rstan","MASS","nlme","ggplot2","phylolm","latex2exp","RColorBrewer","xtable","ggforestplot", "basicPlotteR")

for(lib in LIBRARIES){
  #check if package is installed, if not install and load
  if(!require(lib, character.only = TRUE)){
    install.packages(lib)
    library(lib, character.only = TRUE)
  }
}

if(!require("basicPlotteR")){
  remotes::install_github("JosephCrispell/basicPlotteR")
  library("basicPlotteR")
}

if(!require("ggforestplot")){
  devtools::install_github("NightingaleHealth/ggforestplot")
  library("ggforestplot")
}

## (this just globally sets the ggplot2 theme to theme_bw) 
ggplot2::theme_set(ggplot2::theme_bw(base_size = 16)) 

source("load_data.R")
source("support_functions.R")
source("construct_lifetables.R")
source("primate_tree_bootstap.R");
source("plot_survival_hazard_estimates.R")
source("lifetable_regression_brms.R")
source("gompertzFitLogLikUnion.R");
source("aging_parameters_forest_plot.R")
source("phylo_regression_ancestal_states.R")
source("bayesian_body_weight.R")
source("power_sim2.R")