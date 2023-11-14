
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Define the spatial prior
spatial <- "ICAR"
# spatial <- "PCAR"
# spatial <- "BYM2"



####################
# Load the results #
####################

if(spatial=="ICAR") {
  load("results/MICAR_dowry_rapes_vs_sexratio_2011.Rdata")
  MODELS <- list(micar, micar.eigen64, micar.eigen59, micar.eigen54, micar.eigen49)
}

if(spatial=="PCAR") {
  load("results/MPCAR_dowry_rapes_vs_sexratio_2011.Rdata")
  MODELS <- list(mpcar, mpcar.eigen64, mpcar.eigen59, mpcar.eigen54, mpcar.eigen49)
}

if(spatial=="BYM2") {
  load("results/MBYM2_dowry_rapes_vs_sexratio_2011.Rdata")
  MODELS <- list(mbym2, mbym2.eigen64, mbym2.eigen59, mbym2.eigen54, mbym2.eigen49)
}

names(MODELS) <- c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49")



#############################################################################
# Table 2: posterior mean, sd and 95% credible intervals of beta1 and beta2 #
#############################################################################
Tab.beta1 <- do.call(rbind, lapply(MODELS, function(x) round(x$summary.fixed[3, c(1:3,5)], 4)))
rownames(Tab.beta1) <- names(MODELS)
print(Tab.beta1)


Tab.beta2 <- do.call(rbind, lapply(MODELS, function(x) round(x$summary.fixed[4, c(1:3,5)], 4)))
rownames(Tab.beta2) <- names(MODELS)
print(Tab.beta2)



##########################################################################################
# Table 3: posterior median and 95% credible intervals of the correlation between crimes #
##########################################################################################
cor.crimes <- function(x){
  data.frame(median=x$summary.cor[4], 
             ci.low=x$summary.cor[3],
             ci.up=x$summary.cor[5])}  


Tab.cor <- do.call(rbind, lapply(MODELS, function(x) cor.crimes(x)))
Tab.cor[, 1:3] <- round(Tab.cor[, 1:3], 4)
rownames(Tab.cor) <- names(MODELS)
colnames(Tab.cor) <- c("median", "CI.L", "CI.U")
print(Tab.cor)



#########################################
# Table 4: DIC and WAIC of the M-models #
#########################################
DIC.WAIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance, ## posterior mean deviance
             pD=x$dic$p.eff,                    ## effective number of parameters
             DIC=x$dic$dic,                     ## Deviance Information Criterion
             WAIC=x$waic$waic)}                 ## Watanabe-Akaike information criterion


Tab.DIC.WAIC <- do.call(rbind, lapply(MODELS, DIC.WAIC))
Tab.DIC.WAIC <- round(Tab.DIC.WAIC,4)
rownames(Tab.DIC.WAIC) <- names(MODELS)
print(Tab.DIC.WAIC)




