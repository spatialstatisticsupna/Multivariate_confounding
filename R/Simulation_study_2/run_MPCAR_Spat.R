
################################################################################
#                           Spatial M-models: MPCAR                            #
################################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(INLA)
inla.setOption(inla.mode="classic")
library(spdep)
library(spatialreg)
library(fastDummies)
library(sf)


n.sim <- 300

# Define the scenario
scenario <- "Scenario1"
# scenario <- "Scenario2"
# scenario <- "Scenario3"
# scenario <- "Scenario4"
# scenario <- "Scenario5"



#######################
# Load simulated data #
#######################
load(paste0("Simulated_data/SimuStudy2_", scenario, ".Rdata"))


J <- length(unique(data$Crime))
S <- length(unique(data$dist))



##########################
# Load M-model functions #
##########################
source("functions/mpcar_bartlett.R") 
source("functions/MCAR_single_covariate.R") 



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
if(!file.exists(paste0("results_", scenario)))
      dir.create(paste0("results_", scenario))


for (i in 1:n.sim){
  print(paste0("i=", i))
  micar.spat <- NULL
  # Simulated counts
  data$obs <- simu.O[[i]]
  # Fit the model
  mpcar.spat <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                               O="obs", E="exp", covariate="X1.aux", spatial="proper.bartlett",
                              strategy="simplified.laplace")
  
  data$obs <- NULL

  save(mpcar.spat, file=paste0("results_", scenario, "/MPCAR_Spat_",i,".Rdata"))
}



