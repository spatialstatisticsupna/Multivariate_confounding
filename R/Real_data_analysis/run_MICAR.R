
################################################################################
#                           Spatial M-models: MICAR                            #
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



#################
# Load the data #
#################
load("../../Data/data_UttarPradesh_2011.RData")

J <- length(unique(data$Crime))
S <- length(unique(data$dist))



###########################
# Load M-models functions #
###########################
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


coef <- solve(eigen.vect, data$X1[1:S])


data$X1.eigen64 <- as.vector(scale(eigen.vect[, 1:64]%*%coef[1:64]))
data$X1.eigen59 <- as.vector(scale(eigen.vect[, 1:59]%*%coef[1:59]))
data$X1.eigen54 <- as.vector(scale(eigen.vect[, 1:54]%*%coef[1:54]))
data$X1.eigen49 <- as.vector(scale(eigen.vect[, 1:49]%*%coef[1:49]))



############################
# Multivariate ICAR models #
############################
micar <- MCAR.Covariate(carto=carto_UP, data=data, ID.area="dist", ID.disease="Crime",
                        O="obs", E="exp", covariate="X1", spatial="intrinsic.bartlett",
                        strategy="simplified.laplace")

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


# Folder to save results
if(!file.exists("results")) dir.create("results")

# Save the results
save(micar, micar.eigen64, micar.eigen59, micar.eigen54, micar.eigen49,
     file="results/MICAR_dowry_rapes_vs_sexratio_2011.Rdata")



