
################################################################################
#                           Spatial M-models: MICAR                            #
################################################################################

rm(list=ls())
setwd("")

# Load packages
library(INLA)
inla.setOption(inla.mode="classic")
library(spdep)
library(spatialreg)
library(fastDummies)
library(sf)


n.sim <- 300

# Define the scenario
Scenario <- 1



#######################
# Load simulated data #
#######################
load(paste0("Simulated_data/SimuStudy2_Scenario", Scenario, ".Rdata"))


J <- length(unique(data$Crime))
S <- length(unique(data$dist))



##########################
# Load M-model functions #
##########################
source("R/Simulation_study_2/functions/micar_bartlett.R") 
source("R/Simulation_study_2/functions/MCAR_single_covariate.R") 



###################################
# Define spatial structure matrix #
###################################
nb.mat <- nb2mat(poly2nb(carto_UP), style ="B")
nbs <- unlist(lapply(poly2nb(carto_UP), function(x) length(x)))

Q.xi <- diag(nbs)-nb.mat



#######################################
# Multivariate ICAR models: M-Spatial #
#######################################

# Folder to save results
if(!file.exists(paste0("R/Simulation_study_2/results_Scenario", Scenario)))
      dir.create(paste0("R/Simulation_study_2/results_Scenario", Scenario))


for (i in 1:n.sim){
  print(paste0("i=", i))
  micar.spat <- NULL
  # Simulated counts
  data$obs <- simu.O[[i]]
  # Fit the model
  micar.spat <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                               O="obs", E="exp", covariate="X1.aux", spatial="intrinsic.bartlett",
                              strategy="simplified.laplace")
  
  data$obs <- NULL

  save(micar.spat, file=paste0("R/Simulation_study_2/results_Scenario", Scenario, "/MICAR_Spat_",i,".Rdata"))
}



