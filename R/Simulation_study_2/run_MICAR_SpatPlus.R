
################################################################################
#                          Spatial+ M-models: MICAR                            #
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
source("functions/micar_bartlett.R") 
source("functions/MCAR_single_covariate.R") 



###################################
# Define spatial structure matrix #
###################################
nb.mat <- nb2mat(poly2nb(carto_UP), style ="B")
nbs <- unlist(lapply(poly2nb(carto_UP), function(x) length(x)))

Q.xi <- diag(nbs)-nb.mat



##############################################################
# Sex ratio as a linear combination of the eigenvectors of Q #
##############################################################

# Eigen decomposition of the adjacency matrix
eigendecomp <- eigen(Q.xi)
eigendecomp$values
eigen.vect <- eigendecomp$vectors


coef <- solve(eigen.vect, data$X1.aux[1:S])


data$X1.eigen64 <- as.vector(scale(eigen.vect[, 1:64]%*%coef[1:64]))
data$X1.eigen59 <- as.vector(scale(eigen.vect[, 1:59]%*%coef[1:59]))
data$X1.eigen54 <- as.vector(scale(eigen.vect[, 1:54]%*%coef[1:54]))
data$X1.eigen49 <- as.vector(scale(eigen.vect[, 1:49]%*%coef[1:49]))
data$X1.eigen44 <- as.vector(scale(eigen.vect[, 1:44]%*%coef[1:44]))
data$X1.eigen39 <- as.vector(scale(eigen.vect[, 1:39]%*%coef[1:39]))



########################################
# Multivariate ICAR models: M-SpatPlus #
########################################

# Folder to save results
if(!file.exists(paste0("results_", scenario)))
  dir.create(paste0("results_", scenario))


for (i in 1:n.sim){
  print(paste0("i=", i))
  
  micar.eigen64 <- NULL
  micar.eigen59 <- NULL
  micar.eigen54 <- NULL
  micar.eigen49 <- NULL
  micar.eigen44 <- NULL
  micar.eigen39 <- NULL

  # Simulated counts
  data$obs <- simu.O[[i]]
  # Fit the model
  micar.eigen64 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen64", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  micar.eigen59 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen59", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  micar.eigen54 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen54", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  micar.eigen49 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen49", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  micar.eigen44 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen44", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  micar.eigen39 <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                                  O="obs", E="exp", covariate="X1.eigen39", spatial="intrinsic.bartlett",
                                  strategy="simplified.laplace")
  
  data$obs <- NULL
  
  save(micar.eigen64, micar.eigen59, micar.eigen54, micar.eigen49, micar.eigen44, micar.eigen39,
       file=paste0("results_", scenario, "/MICAR_SpatPlus_",i,".Rdata"))
}





