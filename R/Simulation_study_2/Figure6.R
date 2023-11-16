
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)
library(ggplot2)
library(latex2exp)


n.sim <- 300 


beta1 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[3,1])
}



##############
# Scenario 1 #
##############

# MICAR
load("results_Scenario1/MICAR_SexRatio_2df.Rdata")
load("results_Scenario1/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta1 <- do.call(rbind, lapply(micar.bart, beta1))[, 1]
micar.eigen64.beta1 <- do.call(rbind, lapply(micar.bart.eigen64, beta1))[, 1]
micar.eigen59.beta1 <- do.call(rbind, lapply(micar.bart.eigen59, beta1))[, 1]
micar.eigen54.beta1 <- do.call(rbind, lapply(micar.bart.eigen54, beta1))[, 1]
micar.eigen49.beta1 <- do.call(rbind, lapply(micar.bart.eigen49, beta1))[, 1]
micar.eigen44.beta1 <- do.call(rbind, lapply(micar.bart.eigen44, beta1))[, 1]
micar.eigen39.beta1 <- do.call(rbind, lapply(micar.bart.eigen39, beta1))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, 
   micar.bart.eigen49, micar.bart.eigen44, micar.bart.eigen39)


# MPCAR
load("results_Scenario1/MPCAR_SexRatio_2df.Rdata")
load("results_Scenario1/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta1 <- do.call(rbind, lapply(mpcar.bart, beta1))[, 1]
mpcar.eigen64.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta1))[, 1]
mpcar.eigen59.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta1))[, 1]
mpcar.eigen54.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta1))[, 1]
mpcar.eigen49.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta1))[, 1]
mpcar.eigen44.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen44, beta1))[, 1]
mpcar.eigen39.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen39, beta1))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, 
   mpcar.bart.eigen49, mpcar.bart.eigen44, mpcar.bart.eigen39)


# MBYM2
load("results_Scenario1/MBYM2_SexRatio_2df.Rdata")
load("results_Scenario1/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta1 <- do.call(rbind, lapply(mbym2.bart, beta1))[, 1]
mbym2.eigen64.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta1))[, 1]
mbym2.eigen59.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta1))[, 1]
mbym2.eigen54.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta1))[, 1]
mbym2.eigen49.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta1))[, 1]
mbym2.eigen44.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen44, beta1))[, 1]
mbym2.eigen39.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen39, beta1))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, 
   mbym2.bart.eigen49, mbym2.bart.eigen44, mbym2.bart.eigen39)



###################################
# Scenario 1: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim), rep("M-SpatPlus44", n.sim),
           rep("M-SpatPlus39", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 2100), model=model, 
                         beta=c(micar.beta1, micar.eigen64.beta1, micar.eigen59.beta1,
                                micar.eigen54.beta1, micar.eigen49.beta1, micar.eigen44.beta1,
                                micar.eigen39.beta1))

data.mpcar <- data.frame(prior=rep("MPCAR", 2100), model=model, 
                         beta=c(mpcar.beta1, mpcar.eigen64.beta1, mpcar.eigen59.beta1,
                                mpcar.eigen54.beta1, mpcar.eigen49.beta1, mpcar.eigen44.beta1,
                                mpcar.eigen39.beta1))

data.mbym2 <- data.frame(prior=rep("MBYM2", 2100), model=model, 
                         beta=c(mbym2.beta1, mbym2.eigen64.beta1, mbym2.eigen59.beta1,
                                mbym2.eigen54.beta1, mbym2.eigen49.beta1, mbym2.eigen44.beta1,
                                mbym2.eigen39.beta1))



data.Scenario1 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario1$model <- factor(data.Scenario1$model, levels=c("M-Spatial", "M-SpatPlus64", 
                               "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", 
                               "M-SpatPlus39"))
data.Scenario1$prior <- factor(data.Scenario1$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario1$line <- c(rep(0.15, nrow(data.Scenario1)))
data.Scenario1$scenario <- rep("Scenario 1", nrow(data.Scenario1))



##############
# Scenario 2 #
##############

# MICAR
load("results_Scenario2/MICAR_SexRatio_2df.Rdata")
load("results_Scenario2/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta1 <- do.call(rbind, lapply(micar.bart, beta1))[, 1]
micar.eigen64.beta1 <- do.call(rbind, lapply(micar.bart.eigen64, beta1))[, 1]
micar.eigen59.beta1 <- do.call(rbind, lapply(micar.bart.eigen59, beta1))[, 1]
micar.eigen54.beta1 <- do.call(rbind, lapply(micar.bart.eigen54, beta1))[, 1]
micar.eigen49.beta1 <- do.call(rbind, lapply(micar.bart.eigen49, beta1))[, 1]
micar.eigen44.beta1 <- do.call(rbind, lapply(micar.bart.eigen44, beta1))[, 1]
micar.eigen39.beta1 <- do.call(rbind, lapply(micar.bart.eigen39, beta1))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, 
   micar.bart.eigen49, micar.bart.eigen44, micar.bart.eigen39)


# MPCAR
load("results_Scenario2/MPCAR_SexRatio_2df.Rdata")
load("results_Scenario2/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta1 <- do.call(rbind, lapply(mpcar.bart, beta1))[, 1]
mpcar.eigen64.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta1))[, 1]
mpcar.eigen59.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta1))[, 1]
mpcar.eigen54.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta1))[, 1]
mpcar.eigen49.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta1))[, 1]
mpcar.eigen44.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen44, beta1))[, 1]
mpcar.eigen39.beta1 <- do.call(rbind, lapply(mpcar.bart.eigen39, beta1))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, 
   mpcar.bart.eigen49, mpcar.bart.eigen44, mpcar.bart.eigen39)


# MBYM2
load("results_Scenario2/MBYM2_SexRatio_2df.Rdata")
load("results_Scenario2/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta1 <- do.call(rbind, lapply(mbym2.bart, beta1))[, 1]
mbym2.eigen64.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta1))[, 1]
mbym2.eigen59.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta1))[, 1]
mbym2.eigen54.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta1))[, 1]
mbym2.eigen49.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta1))[, 1]
mbym2.eigen44.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen44, beta1))[, 1]
mbym2.eigen39.beta1 <- do.call(rbind, lapply(mbym2.bart.eigen39, beta1))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, 
   mbym2.bart.eigen49, mbym2.bart.eigen44, mbym2.bart.eigen39)



###################################
# Scenario 2: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim), rep("M-SpatPlus44", n.sim),
           rep("M-SpatPlus39", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 2100), model=model, 
                         beta=c(micar.beta1, micar.eigen64.beta1, micar.eigen59.beta1,
                                micar.eigen54.beta1, micar.eigen49.beta1, micar.eigen44.beta1,
                                micar.eigen39.beta1))

data.mpcar <- data.frame(prior=rep("MPCAR", 2100), model=model, 
                         beta=c(mpcar.beta1, mpcar.eigen64.beta1, mpcar.eigen59.beta1,
                                mpcar.eigen54.beta1, mpcar.eigen49.beta1, mpcar.eigen44.beta1,
                                mpcar.eigen39.beta1))

data.mbym2 <- data.frame(prior=rep("MBYM2", 2100), model=model, 
                         beta=c(mbym2.beta1, mbym2.eigen64.beta1, mbym2.eigen59.beta1,
                                mbym2.eigen54.beta1, mbym2.eigen49.beta1, mbym2.eigen44.beta1,
                                mbym2.eigen39.beta1))



data.Scenario2 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario2$model <- factor(data.Scenario2$model, levels=c("M-Spatial", "M-SpatPlus64", 
                               "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", 
                               "M-SpatPlus39"))
data.Scenario2$prior <- factor(data.Scenario2$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario2$line <- c(rep(0.15, nrow(data.Scenario2)))
data.Scenario2$scenario <- rep("Scenario 2", nrow(data.Scenario2))



#########################################
# Merge datasets of different scenarios #
#########################################
data <- rbind(data.Scenario1, data.Scenario2)
data$model <- factor(data$model, levels=c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", 
                    "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", "M-SpatPlus39"))
data$prior <- factor(data$prior, levels=c("MICAR", "MPCAR", "MBYM2"))



######################################################
# Figure 4: boxplot of the estimated means of beta 1 #
######################################################
pdf(file="Figure6.pdf", width=14, height=8)

p <- ggplot(data, aes(x = model, y = beta, fill=prior)) + 
  geom_boxplot(position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9),
        text = element_text(size = 20),
        panel.spacing = unit(1.3, "lines")) + 
  xlab("") + ylab(TeX('$\\beta_{1}$'))  +
  scale_y_continuous(limits=c(0.1, 0.4), breaks=c(0.1,0.2,0.3,0.4))

p + geom_hline(aes(yintercept = c(rep(0.15,nrow(data)))), color="#e31a1c", linetype="dashed") + 
  facet_wrap(~scenario, nrow=2, strip.position="right") + 
  theme(strip.background =element_rect(fill="lightskyblue")) + 
  scale_fill_manual(values = c("#fbb4ae","#b3cde3","#ccebc5")) 


dev.off()

