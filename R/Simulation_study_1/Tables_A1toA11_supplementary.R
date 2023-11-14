
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)


n.sim <- 300 


# Define the scenario
scenario <- "Scenario3"
# scenario <- "Scenario4"
# scenario <- "Scenario5"
# scenario <- "Scenario6"


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



############################################
# Table A1: posterior mean and sd of beta1 #
############################################
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



############################################
# Table A2: posterior mean and sd of beta2 #
############################################
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



###############################################################
# Table A3: estimated and simulated standard errors for beta1 #
###############################################################

# Auxiliar functions
beta1.mean <- function(x){
  data.frame(alpha.mean=x$summary.fixed[3,1])
}

# Sim and est functions
est.beta1 <- function(Model){
  data.frame(apply(do.call(rbind, lapply(Model, beta1))[2], 2, mean))
}

sim.beta1 <- function(Model){
  all.mean <- rep(rbind(apply(do.call(rbind, lapply(Model, beta1))[1], 2, mean)), n.sim)
  data.frame(sqrt(Reduce("+",mapply(function(x,y){(x-y)^2}, 
      x=lapply(Model, beta1.mean), y=all.mean, SIMPLIFY=FALSE))/n.sim))
}


Tab.est.beta1 <- data.frame(Spat=est.beta1(MSpat),
                            SpatPlus64=est.beta1(MSpatPlus64),
                            SpatPlus59=est.beta1(MSpatPlus59),
                            SpatPlus54=est.beta1(MSpatPlus54),
                            SpatPlus49=est.beta1(MSpatPlus49))
colnames(Tab.est.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")


Tab.sim.beta1 <- data.frame(Spat=sim.beta1(MSpat),
                            SpatPlus64=sim.beta1(MSpatPlus64),
                            SpatPlus59=sim.beta1(MSpatPlus59),
                            SpatPlus54=sim.beta1(MSpatPlus54),
                            SpatPlus49=sim.beta1(MSpatPlus49))
colnames(Tab.sim.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")


Tab.sim.est.beta1 <- rbind(Tab.est.beta1, Tab.sim.beta1)
rownames(Tab.sim.est.beta1) <- c("estimated", "simulated")

Tab.sim.est.beta1[, 1:5] <- round(Tab.sim.est.beta1[, 1:5], 4)
Tab.sim.est.beta1 <- data.frame(t(Tab.sim.est.beta1))
print(Tab.sim.est.beta1)



###############################################################
# Table A4: estimated and simulated standard errors for beta2 #
###############################################################

# Auxiliar functions
beta2.mean <- function(x){
  data.frame(alpha.mean=x$summary.fixed[4,1])
}

# Sim and est functions
est.beta2 <- function(Model){
  data.frame(apply(do.call(rbind, lapply(Model, beta2))[2], 2, mean))
}

sim.beta2 <- function(Model){
  all.mean <- rep(rbind(apply(do.call(rbind, lapply(Model, beta2))[1], 2, mean)), n.sim)
  data.frame(sqrt(Reduce("+",mapply(function(x,y){(x-y)^2}, 
                                    x=lapply(Model, beta2.mean), y=all.mean, SIMPLIFY=FALSE))/n.sim))
}


Tab.est.beta2 <- data.frame(Spat=est.beta2(MSpat),
                            SpatPlus64=est.beta2(MSpatPlus64),
                            SpatPlus59=est.beta2(MSpatPlus59),
                            SpatPlus54=est.beta2(MSpatPlus54),
                            SpatPlus49=est.beta2(MSpatPlus49))
colnames(Tab.est.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")


Tab.sim.beta2 <- data.frame(Spat=sim.beta2(MSpat),
                            SpatPlus64=sim.beta2(MSpatPlus64),
                            SpatPlus59=sim.beta2(MSpatPlus59),
                            SpatPlus54=sim.beta2(MSpatPlus54),
                            SpatPlus49=sim.beta2(MSpatPlus49))
colnames(Tab.sim.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")


Tab.sim.est.beta2 <- rbind(Tab.est.beta2, Tab.sim.beta2)
rownames(Tab.sim.est.beta2) <- c("estimated", "simulated")

Tab.sim.est.beta2[, 1:5] <- round(Tab.sim.est.beta2[, 1:5], 4)
Tab.sim.est.beta2 <- data.frame(t(Tab.sim.est.beta2))
print(Tab.sim.est.beta2)



####################################################################
# Table A5a: 95% coverage probabilities of the true value of beta1 #
####################################################################
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



####################################################################
# Table A5b: 95% coverage probabilities of the true value of beta2 #
####################################################################
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



###########################################################################
# Table A6: posterior median and 95% CI of the correlation between crimes #
###########################################################################
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



########################################
# Table A7: DIC and WAIC of the models #
########################################
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



###################################################
# Table A8: MARB and MRRMSE of the relative risks #
###################################################
MARB.RR <- function(Model){
  data.frame(MARB=mean(abs(Reduce("+",mapply(function(x,y){(x-y)/y}, 
             x=Model, y=Real_risk_vect, SIMPLIFY=FALSE)))/n.sim))
}

MRRMSE.RR <- function(Model){
  data.frame(mean(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, 
             x=Model, y=Real_risk_vect, SIMPLIFY=FALSE))/n.sim)))
}


# Load the simulated data
load(paste0("Simulated_data/SimuStudy1_", scenario, ".Rdata"))

Real.risk <- exp(log.risk)

Real_risk_vect <- vector("list", n.sim)
for (i in 1:n.sim){
  Real_risk_vect[[i]] <- Real.risk
}


RR.Spat <- lapply(MSpat, function(x) matrix(x$summary.fitted.values$'0.5quant',140,1))
RR.SpatPlus64 <- lapply(MSpatPlus64, function(x) matrix(x$summary.fitted.values$'0.5quant',140,1))
RR.SpatPlus59 <- lapply(MSpatPlus59, function(x) matrix(x$summary.fitted.values$'0.5quant',140,1))
RR.SpatPlus54 <- lapply(MSpatPlus54, function(x) matrix(x$summary.fitted.values$'0.5quant',140,1))
RR.SpatPlus49 <- lapply(MSpatPlus49, function(x) matrix(x$summary.fitted.values$'0.5quant',140,1))



MARB <- data.frame(Spat=MARB.RR(RR.Spat),
                   SpatPlus64=MARB.RR(RR.SpatPlus64),
                   SpatPlus59=MARB.RR(RR.SpatPlus59),
                   SpatPlus54=MARB.RR(RR.SpatPlus54),
                   SpatPlus49=MARB.RR(RR.SpatPlus49))


colnames(MARB) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MARB[, 1:5] <- round(MARB[, 1:5], 4)
MARB <- data.frame(t(MARB))
colnames(MARB) <- "MARB"
print(MARB)


MRRMSE <- data.frame(Spat=MRRMSE.RR(RR.Spat),
                     SpatPlus64=MRRMSE.RR(RR.SpatPlus64),
                     SpatPlus59=MRRMSE.RR(RR.SpatPlus59),
                     SpatPlus54=MRRMSE.RR(RR.SpatPlus54),
                     SpatPlus49=MRRMSE.RR(RR.SpatPlus49))


colnames(MRRMSE) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MRRMSE[, 1:5] <- round(MRRMSE[, 1:5], 4)
MRRMSE <- data.frame(t(MRRMSE))
colnames(MRRMSE) <- "MRRMSE"
print(MRRMSE)



Table.Accuracy <- cbind(MARB, MRRMSE)
Table.Accuracy[, 1:2] <- round(Table.Accuracy[, 1:2], 4)
Table.Accuracy <- data.frame(Table.Accuracy)
print(Table.Accuracy)



#################################################################
# Table A9: length of 95% credible intervals of beta1 and beta2 # 
#################################################################
Length.beta1 <- function(Model){
  data.frame(length=abs(Model$summary.fixed$'0.975quant'[3]-Model$summary.fixed$'0.025quant'[3]))
}

Length.beta2 <- function(Model){
  data.frame(length=abs(Model$summary.fixed$'0.975quant'[4]-Model$summary.fixed$'0.025quant'[4]))
}



Length.beta1 <- data.frame(Spat=mean(unlist(lapply(MSpat, function(x) Length.beta1(x)))),
                           SpatPlus64=mean(unlist(lapply(MSpatPlus64, function(x) Length.beta1(x)))),
                           SpatPlus59=mean(unlist(lapply(MSpatPlus59, function(x) Length.beta1(x)))),
                           SpatPlus54=mean(unlist(lapply(MSpatPlus54, function(x) Length.beta1(x)))),
                           SpatPlus49=mean(unlist(lapply(MSpatPlus49, function(x) Length.beta1(x)))))   

Length.beta1 <- as.data.frame(t(Length.beta1))
rownames(Length.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Length.beta1) <- c("beta1")
print(Length.beta1)


Length.beta2 <- data.frame(Spat=mean(unlist(lapply(MSpat, function(x) Length.beta2(x)))),
                           SpatPlus64=mean(unlist(lapply(MSpatPlus64, function(x) Length.beta2(x)))),
                           SpatPlus59=mean(unlist(lapply(MSpatPlus59, function(x) Length.beta2(x)))),
                           SpatPlus54=mean(unlist(lapply(MSpatPlus54, function(x) Length.beta2(x)))),
                           SpatPlus49=mean(unlist(lapply(MSpatPlus49, function(x) Length.beta2(x)))))   

Length.beta2 <- as.data.frame(t(Length.beta2))
rownames(Length.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
colnames(Length.beta2) <- c("beta2")
print(Length.beta2)



#######################################
# Table A10: MARB and MRRMSE of beta1 #
#######################################
real.beta1 <- -0.15
real.beta1.vect <- rep(real.beta1, n.sim)


beta1.Spat <- lapply(MSpat, beta1.mean)
beta1.SpatPlus64 <- lapply(MSpatPlus64, beta1.mean)
beta1.SpatPlus59 <- lapply(MSpatPlus59, beta1.mean)
beta1.SpatPlus54 <- lapply(MSpatPlus54, beta1.mean)
beta1.SpatPlus49 <- lapply(MSpatPlus49, beta1.mean)


MARB.beta1 <- function(Model){
  data.frame(MARB=abs(Reduce("+",mapply(function(x,y){(x-y)/y}, 
                                        x=Model, y=real.beta1.vect, SIMPLIFY=FALSE)))/n.sim)
}


MRRMSE.beta1 <- function(Model){
  data.frame(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, 
                                    x=Model, y=real.beta1.vect, SIMPLIFY=FALSE))/n.sim))
}


MARB.beta1 <- data.frame(Spat=MARB.beta1(beta1.Spat),
                         SpatPlus64=MARB.beta1(beta1.SpatPlus64),
                         SpatPlus59=MARB.beta1(beta1.SpatPlus59),
                         SpatPlus54=MARB.beta1(beta1.SpatPlus54),
                         SpatPlus49=MARB.beta1(beta1.SpatPlus49))


colnames(MARB.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MARB.beta1[, 1:5] <- round(MARB.beta1[, 1:5], 4)
MARB.beta1 <- data.frame(t(MARB.beta1))
colnames(MARB.beta1) <- "MARB.beta1"
print(MARB.beta1)



MRRMSE.beta1 <- data.frame(Spat=MRRMSE.beta1(beta1.Spat),
                           SpatPlus64=MRRMSE.beta1(beta1.SpatPlus64),
                           SpatPlus59=MRRMSE.beta1(beta1.SpatPlus59),
                           SpatPlus54=MRRMSE.beta1(beta1.SpatPlus54),
                           SpatPlus49=MRRMSE.beta1(beta1.SpatPlus49))


colnames(MRRMSE.beta1) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MRRMSE.beta1[, 1:5] <- round(MRRMSE.beta1[, 1:5], 4)
MRRMSE.beta1 <- data.frame(t(MRRMSE.beta1))
colnames(MRRMSE.beta1) <- "MRRMSE.beta1"
print(MRRMSE.beta1)



#######################################
# Table A11: MARB and MRRMSE of beta2 #
#######################################
real.beta2 <- -0.2
real.beta2.vect <- rep(real.beta2, n.sim)


beta2.Spat <- lapply(MSpat, beta2.mean)
beta2.SpatPlus64 <- lapply(MSpatPlus64, beta2.mean)
beta2.SpatPlus59 <- lapply(MSpatPlus59, beta2.mean)
beta2.SpatPlus54 <- lapply(MSpatPlus54, beta2.mean)
beta2.SpatPlus49 <- lapply(MSpatPlus49, beta2.mean)


MARB.beta2 <- function(Model){
  data.frame(MARB=abs(Reduce("+",mapply(function(x,y){(x-y)/y}, 
                                        x=Model, y=real.beta2.vect, SIMPLIFY=FALSE)))/n.sim)
}


MRRMSE.beta2 <- function(Model){
  data.frame(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, 
                                    x=Model, y=real.beta2.vect, SIMPLIFY=FALSE))/n.sim))
}


MARB.beta2 <- data.frame(Spat=MARB.beta2(beta2.Spat),
                         SpatPlus64=MARB.beta2(beta2.SpatPlus64),
                         SpatPlus59=MARB.beta2(beta2.SpatPlus59),
                         SpatPlus54=MARB.beta2(beta2.SpatPlus54),
                         SpatPlus49=MARB.beta2(beta2.SpatPlus49))


colnames(MARB.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MARB.beta2[, 1:5] <- round(MARB.beta2[, 1:5], 4)
MARB.beta2 <- data.frame(t(MARB.beta2))
colnames(MARB.beta2) <- "MARB.beta2"
print(MARB.beta2)



MRRMSE.beta2 <- data.frame(Spat=MRRMSE.beta2(beta2.Spat),
                           SpatPlus64=MRRMSE.beta2(beta2.SpatPlus64),
                           SpatPlus59=MRRMSE.beta2(beta2.SpatPlus59),
                           SpatPlus54=MRRMSE.beta2(beta2.SpatPlus54),
                           SpatPlus49=MRRMSE.beta2(beta2.SpatPlus49))


colnames(MRRMSE.beta2) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")
MRRMSE.beta2[, 1:5] <- round(MRRMSE.beta2[, 1:5], 4)
MRRMSE.beta2 <- data.frame(t(MRRMSE.beta2))
colnames(MRRMSE.beta2) <- "MRRMSE.beta2"
print(MRRMSE.beta2)


