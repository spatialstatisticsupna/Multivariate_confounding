
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)


n.sim <- 300 


# Define the scenario
scenario <- "Scenario1"
# scenario <- "Scenario2"


# Choose the spatial prior
spatial <- "MICAR"
# spatial <- "MPCAR"
# spatial <- "MBYM2"



####################
# Load the results #
####################
load(paste0("results_", scenario, "/", spatial,"_SexRatio_2df.Rdata"))
load(paste0("results_", scenario, "/", spatial,"_SpatPlus_SexRatio_2df.Rdata"))


if (spatial=="MICAR"){
  MSpat <- micar.bart
  MSpatPlus64 <- micar.bart.eigen64
  MSpatPlus59 <- micar.bart.eigen59
  MSpatPlus54 <- micar.bart.eigen54
  MSpatPlus49 <- micar.bart.eigen49
}

if (spatial=="MPCAR"){
  MSpat <- mpcar.bart
  MSpatPlus64 <- mpcar.bart.eigen64
  MSpatPlus59 <- mpcar.bart.eigen59
  MSpatPlus54 <- mpcar.bart.eigen54
  MSpatPlus49 <- mpcar.bart.eigen49
}

if (spatial=="MBYM2"){
  MSpat <- mbym2.bart
  MSpatPlus64 <- mbym2.bart.eigen64
  MSpatPlus59 <- mbym2.bart.eigen59
  MSpatPlus54 <- mbym2.bart.eigen54
  MSpatPlus49 <- mbym2.bart.eigen49
}



###########################################
# Table 5: posterior mean and sd of beta1 #
###########################################
beta1 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[3,1], 
             alpha.sd=x$summary.fixed[3,2])
}


Spat.beta1 <- rbind(apply(do.call(rbind, lapply(MSpat, beta1)), 2, mean))
SpatPlus64.beta1 <- rbind(apply(do.call(rbind, lapply(MSpatPlus64, beta1)), 2, mean))
SpatPlus59.beta1 <- rbind(apply(do.call(rbind, lapply(MSpatPlus59, beta1)), 2, mean))
SpatPlus54.beta1 <- rbind(apply(do.call(rbind, lapply(MSpatPlus54, beta1)), 2, mean))
SpatPlus49.beta1 <- rbind(apply(do.call(rbind, lapply(MSpatPlus49, beta1)), 2, mean))


Tab.beta1 <- rbind(Spat.beta1, SpatPlus64.beta1, SpatPlus59.beta1, SpatPlus54.beta1, SpatPlus49.beta1)
Tab.beta1 <- as.data.frame(Tab.beta1)
rownames(Tab.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Tab.beta1) <- c("mean", "sd")
Tab.beta1[, 1:2] <- round(Tab.beta1[, 1:2], 4)
print(Tab.beta1)



###########################################
# Table 6: posterior mean and sd of beta2 #
###########################################
beta2 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[4,1], 
             alpha.sd=x$summary.fixed[4,2])
}


Spat.beta2 <- rbind(apply(do.call(rbind, lapply(MSpat, beta2)), 2, mean))
SpatPlus64.beta2 <- rbind(apply(do.call(rbind, lapply(MSpatPlus64, beta2)), 2, mean))
SpatPlus59.beta2 <- rbind(apply(do.call(rbind, lapply(MSpatPlus59, beta2)), 2, mean))
SpatPlus54.beta2 <- rbind(apply(do.call(rbind, lapply(MSpatPlus54, beta2)), 2, mean))
SpatPlus49.beta2 <- rbind(apply(do.call(rbind, lapply(MSpatPlus49, beta2)), 2, mean))


Tab.beta2 <- rbind(Spat.beta2, SpatPlus64.beta2, SpatPlus59.beta2, SpatPlus54.beta2, SpatPlus49.beta2)
Tab.beta2 <- as.data.frame(Tab.beta2)
rownames(Tab.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Tab.beta2) <- c("mean", "sd")
Tab.beta2[, 1:2] <- round(Tab.beta2[, 1:2], 4)
print(Tab.beta2)



###################################################################
# Table 7a: 95% coverage probabilities of the true value of beta1 #
###################################################################
real.beta1 <- -0.15

Coverage.beta1 <- function(Model, true.value){
  CI.L <- Model$summary.fixed$'0.025quant'[3]
  CI.U <- Model$summary.fixed$'0.975quant'[3]
  ifelse(true.value>=CI.L & true.value<=CI.U, 1, 0)
}


Coverage.beta1 <- data.frame(Spat=mean(unlist(lapply(MSpat, function(x) Coverage.beta1(x,real.beta1))))*100,
                             SpatPlus64=mean(unlist(lapply(MSpatPlus64, function(x) Coverage.beta1(x,real.beta1))))*100,
                             SpatPlus59=mean(unlist(lapply(MSpatPlus59, function(x) Coverage.beta1(x,real.beta1))))*100,
                             SpatPlus54=mean(unlist(lapply(MSpatPlus54, function(x) Coverage.beta1(x,real.beta1))))*100,
                             SpatPlus49=mean(unlist(lapply(MSpatPlus49, function(x) Coverage.beta1(x,real.beta1))))*100)   


Coverage.beta1 <- as.data.frame(t(Coverage.beta1))
Coverage.beta1[, 1] <- round(Coverage.beta1[, 1], 4)
rownames(Coverage.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Coverage.beta1) <- paste0(spatial, ".beta1")
print(Coverage.beta1)



###################################################################
# Table 7b: 95% coverage probabilities of the true value of beta2 #
###################################################################
real.beta2 <- -0.2

Coverage.beta2 <- function(Model, true.value){
  CI.L <- Model$summary.fixed$'0.025quant'[4]
  CI.U <- Model$summary.fixed$'0.975quant'[4]
  ifelse(true.value>=CI.L & true.value<=CI.U, 1, 0)
}


Coverage.beta2 <- data.frame(Spat=mean(unlist(lapply(MSpat, function(x) Coverage.beta2(x,real.beta2))))*100,
                             SpatPlus64=mean(unlist(lapply(MSpatPlus64, function(x) Coverage.beta2(x,real.beta2))))*100,
                             SpatPlus59=mean(unlist(lapply(MSpatPlus59, function(x) Coverage.beta2(x,real.beta2))))*100,
                             SpatPlus54=mean(unlist(lapply(MSpatPlus54, function(x) Coverage.beta2(x,real.beta2))))*100,
                             SpatPlus49=mean(unlist(lapply(MSpatPlus49, function(x) Coverage.beta2(x,real.beta2))))*100)   


Coverage.beta2 <- as.data.frame(t(Coverage.beta2))
Coverage.beta2[, 1] <- round(Coverage.beta2[, 1], 4)
rownames(Coverage.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Coverage.beta2) <- paste0(spatial, ".beta2")
print(Coverage.beta2)



##########################################################################
# Table 8: posterior median and 95% CI of the correlation between crimes #
##########################################################################
cor.crimes <- function(x){
  data.frame(median=x$summary.cor[4], 
             ci.low=x$summary.cor[3],
             ci.up=x$summary.cor[5])}  


Spat.cor <- rbind(apply(do.call(rbind, lapply(MSpat, cor.crimes)), 2, mean))
SpatPlus64.cor <- rbind(apply(do.call(rbind, lapply(MSpatPlus64, cor.crimes)), 2, mean))
SpatPlus59.cor <- rbind(apply(do.call(rbind, lapply(MSpatPlus59, cor.crimes)), 2, mean))
SpatPlus54.cor <- rbind(apply(do.call(rbind, lapply(MSpatPlus54, cor.crimes)), 2, mean))
SpatPlus49.cor <- rbind(apply(do.call(rbind, lapply(MSpatPlus49, cor.crimes)), 2, mean))


Tab.cor <- rbind(Spat.cor, SpatPlus64.cor, SpatPlus59.cor, SpatPlus54.cor, SpatPlus49.cor)
Tab.cor <- as.data.frame(Tab.cor)
rownames(Tab.cor) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Tab.cor) <- c("median", "CI.L", "CI.U")
Tab.cor[, 1:3] <- round(Tab.cor[, 1:3], 4)
print(Tab.cor)



#######################################
# Table 9: DIC and WAIC of the models #
#######################################
DIC.WAIC <- function(x){
  data.frame(DIC=x$dic$dic, WAIC=x$waic$waic)
}                 


Spat.dic.waic <- rbind(apply(do.call(rbind, lapply(MSpat, DIC.WAIC)), 2, mean))
SpatPlus64.dic.waic <- rbind(apply(do.call(rbind, lapply(MSpatPlus64, DIC.WAIC)), 2, mean))
SpatPlus59.dic.waic <- rbind(apply(do.call(rbind, lapply(MSpatPlus59, DIC.WAIC)), 2, mean))
SpatPlus54.dic.waic <- rbind(apply(do.call(rbind, lapply(MSpatPlus54, DIC.WAIC)), 2, mean))
SpatPlus49.dic.waic <- rbind(apply(do.call(rbind, lapply(MSpatPlus49, DIC.WAIC)), 2, mean))


Tab.DIC.WAIC <- rbind(Spat.dic.waic, SpatPlus64.dic.waic, SpatPlus59.dic.waic, 
                      SpatPlus54.dic.waic, SpatPlus49.dic.waic)
Tab.DIC.WAIC <- as.data.frame(Tab.DIC.WAIC)
rownames(Tab.DIC.WAIC) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Tab.DIC.WAIC) <- c("DIC", "WAIC")
Tab.DIC.WAIC[, 1:2] <- round(Tab.DIC.WAIC[, 1:2], 4)
print(Tab.DIC.WAIC)






