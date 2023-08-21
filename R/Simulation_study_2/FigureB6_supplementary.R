

rm(list=ls())
setwd("")


# Load packages
library(INLA)
library(ggplot2)
library(latex2exp)


n.sim <- 300 


beta2 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[4,1])
}



##############
# Scenario 3 #
##############

# MICAR
load("Simulation_study_2/results_Scenario3/MICAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario3/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta2 <- do.call(rbind, lapply(micar.bart, beta2))[, 1]
micar.eigen64.beta2 <- do.call(rbind, lapply(micar.bart.eigen64, beta2))[, 1]
micar.eigen59.beta2 <- do.call(rbind, lapply(micar.bart.eigen59, beta2))[, 1]
micar.eigen54.beta2 <- do.call(rbind, lapply(micar.bart.eigen54, beta2))[, 1]
micar.eigen49.beta2 <- do.call(rbind, lapply(micar.bart.eigen49, beta2))[, 1]
micar.eigen44.beta2 <- do.call(rbind, lapply(micar.bart.eigen44, beta2))[, 1]
micar.eigen39.beta2 <- do.call(rbind, lapply(micar.bart.eigen39, beta2))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, 
   micar.bart.eigen49, micar.bart.eigen44, micar.bart.eigen39)


# MPCAR
load("Simulation_study_2/results_Scenario3/MPCAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario3/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta2 <- do.call(rbind, lapply(mpcar.bart, beta2))[, 1]
mpcar.eigen64.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta2))[, 1]
mpcar.eigen59.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta2))[, 1]
mpcar.eigen54.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta2))[, 1]
mpcar.eigen49.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta2))[, 1]
mpcar.eigen44.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen44, beta2))[, 1]
mpcar.eigen39.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen39, beta2))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, 
   mpcar.bart.eigen49, mpcar.bart.eigen44, mpcar.bart.eigen39)


# MBYM2
load("Simulation_study_2/results_Scenario3/MBYM2_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario3/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta2 <- do.call(rbind, lapply(mbym2.bart, beta2))[, 1]
mbym2.eigen64.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta2))[, 1]
mbym2.eigen59.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta2))[, 1]
mbym2.eigen54.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta2))[, 1]
mbym2.eigen49.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta2))[, 1]
mbym2.eigen44.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen44, beta2))[, 1]
mbym2.eigen39.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen39, beta2))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, 
   mbym2.bart.eigen49, mbym2.bart.eigen44, mbym2.bart.eigen39)



###################################
# Scenario 3: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim), rep("M-SpatPlus44", n.sim),
           rep("M-SpatPlus39", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 2100), model=model, 
                         beta=c(micar.beta2, micar.eigen64.beta2, micar.eigen59.beta2,
                                micar.eigen54.beta2, micar.eigen49.beta2, micar.eigen44.beta2,
                                micar.eigen39.beta2))

data.mpcar <- data.frame(prior=rep("MPCAR", 2100), model=model, 
                         beta=c(mpcar.beta2, mpcar.eigen64.beta2, mpcar.eigen59.beta2,
                                mpcar.eigen54.beta2, mpcar.eigen49.beta2, mpcar.eigen44.beta2,
                                mpcar.eigen39.beta2))

data.mbym2 <- data.frame(prior=rep("MBYM2", 2100), model=model, 
                         beta=c(mbym2.beta2, mbym2.eigen64.beta2, mbym2.eigen59.beta2,
                                mbym2.eigen54.beta2, mbym2.eigen49.beta2, mbym2.eigen44.beta2,
                                mbym2.eigen39.beta2))



data.Scenario3 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario3$model <- factor(data.Scenario3$model, levels=c("M-Spatial", "M-SpatPlus64", 
                               "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", 
                               "M-SpatPlus39"))
data.Scenario3$prior <- factor(data.Scenario3$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario3$line <- c(rep(0.2, nrow(data.Scenario3)))
data.Scenario3$scenario <- rep("Scenario 3", nrow(data.Scenario3))



##############
# Scenario 4 #
##############

# MICAR
load("Simulation_study_2/results_Scenario4/MICAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario4/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta2 <- do.call(rbind, lapply(micar.bart, beta2))[, 1]
micar.eigen64.beta2 <- do.call(rbind, lapply(micar.bart.eigen64, beta2))[, 1]
micar.eigen59.beta2 <- do.call(rbind, lapply(micar.bart.eigen59, beta2))[, 1]
micar.eigen54.beta2 <- do.call(rbind, lapply(micar.bart.eigen54, beta2))[, 1]
micar.eigen49.beta2 <- do.call(rbind, lapply(micar.bart.eigen49, beta2))[, 1]
micar.eigen44.beta2 <- do.call(rbind, lapply(micar.bart.eigen44, beta2))[, 1]
micar.eigen39.beta2 <- do.call(rbind, lapply(micar.bart.eigen39, beta2))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, 
   micar.bart.eigen49, micar.bart.eigen44, micar.bart.eigen39)


# MPCAR
load("Simulation_study_2/results_Scenario4/MPCAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario4/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta2 <- do.call(rbind, lapply(mpcar.bart, beta2))[, 1]
mpcar.eigen64.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta2))[, 1]
mpcar.eigen59.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta2))[, 1]
mpcar.eigen54.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta2))[, 1]
mpcar.eigen49.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta2))[, 1]
mpcar.eigen44.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen44, beta2))[, 1]
mpcar.eigen39.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen39, beta2))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, 
   mpcar.bart.eigen49, mpcar.bart.eigen44, mpcar.bart.eigen39)


# MBYM2
load("Simulation_study_2/results_Scenario4/MBYM2_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario4/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta2 <- do.call(rbind, lapply(mbym2.bart, beta2))[, 1]
mbym2.eigen64.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta2))[, 1]
mbym2.eigen59.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta2))[, 1]
mbym2.eigen54.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta2))[, 1]
mbym2.eigen49.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta2))[, 1]
mbym2.eigen44.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen44, beta2))[, 1]
mbym2.eigen39.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen39, beta2))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, 
   mbym2.bart.eigen49, mbym2.bart.eigen44, mbym2.bart.eigen39)



###################################
# Scenario 4: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim), rep("M-SpatPlus44", n.sim),
           rep("M-SpatPlus39", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 2100), model=model, 
                         beta=c(micar.beta2, micar.eigen64.beta2, micar.eigen59.beta2,
                                micar.eigen54.beta2, micar.eigen49.beta2, micar.eigen44.beta2,
                                micar.eigen39.beta2))

data.mpcar <- data.frame(prior=rep("MPCAR", 2100), model=model, 
                         beta=c(mpcar.beta2, mpcar.eigen64.beta2, mpcar.eigen59.beta2,
                                mpcar.eigen54.beta2, mpcar.eigen49.beta2, mpcar.eigen44.beta2,
                                mpcar.eigen39.beta2))

data.mbym2 <- data.frame(prior=rep("MBYM2", 2100), model=model, 
                         beta=c(mbym2.beta2, mbym2.eigen64.beta2, mbym2.eigen59.beta2,
                                mbym2.eigen54.beta2, mbym2.eigen49.beta2, mbym2.eigen44.beta2,
                                mbym2.eigen39.beta2))



data.Scenario4 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario4$model <- factor(data.Scenario4$model, levels=c("M-Spatial", "M-SpatPlus64", 
                               "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", 
                               "M-SpatPlus39"))
data.Scenario4$prior <- factor(data.Scenario4$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario4$line <- c(rep(0.2, nrow(data.Scenario4)))
data.Scenario4$scenario <- rep("Scenario 4", nrow(data.Scenario4))



##############
# Scenario 5 #
##############

# MICAR
load("Simulation_study_2/results_Scenario5/MICAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario5/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta2 <- do.call(rbind, lapply(micar.bart, beta2))[, 1]
micar.eigen64.beta2 <- do.call(rbind, lapply(micar.bart.eigen64, beta2))[, 1]
micar.eigen59.beta2 <- do.call(rbind, lapply(micar.bart.eigen59, beta2))[, 1]
micar.eigen54.beta2 <- do.call(rbind, lapply(micar.bart.eigen54, beta2))[, 1]
micar.eigen49.beta2 <- do.call(rbind, lapply(micar.bart.eigen49, beta2))[, 1]
micar.eigen44.beta2 <- do.call(rbind, lapply(micar.bart.eigen44, beta2))[, 1]
micar.eigen39.beta2 <- do.call(rbind, lapply(micar.bart.eigen39, beta2))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, 
   micar.bart.eigen49, micar.bart.eigen44, micar.bart.eigen39)


# MPCAR
load("Simulation_study_2/results_Scenario5/MPCAR_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario5/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta2 <- do.call(rbind, lapply(mpcar.bart, beta2))[, 1]
mpcar.eigen64.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta2))[, 1]
mpcar.eigen59.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta2))[, 1]
mpcar.eigen54.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta2))[, 1]
mpcar.eigen49.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta2))[, 1]
mpcar.eigen44.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen44, beta2))[, 1]
mpcar.eigen39.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen39, beta2))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, 
   mpcar.bart.eigen49, mpcar.bart.eigen44, mpcar.bart.eigen39)


# MBYM2
load("Simulation_study_2/results_Scenario5/MBYM2_SexRatio_2df.Rdata")
load("Simulation_study_2/results_Scenario5/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta2 <- do.call(rbind, lapply(mbym2.bart, beta2))[, 1]
mbym2.eigen64.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta2))[, 1]
mbym2.eigen59.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta2))[, 1]
mbym2.eigen54.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta2))[, 1]
mbym2.eigen49.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta2))[, 1]
mbym2.eigen44.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen44, beta2))[, 1]
mbym2.eigen39.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen39, beta2))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, 
   mbym2.bart.eigen49, mbym2.bart.eigen44, mbym2.bart.eigen39)



###################################
# Scenario 5: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim), rep("M-SpatPlus44", n.sim),
           rep("M-SpatPlus39", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 2100), model=model, 
                         beta=c(micar.beta2, micar.eigen64.beta2, micar.eigen59.beta2,
                                micar.eigen54.beta2, micar.eigen49.beta2, micar.eigen44.beta2,
                                micar.eigen39.beta2))

data.mpcar <- data.frame(prior=rep("MPCAR", 2100), model=model, 
                         beta=c(mpcar.beta2, mpcar.eigen64.beta2, mpcar.eigen59.beta2,
                                mpcar.eigen54.beta2, mpcar.eigen49.beta2, mpcar.eigen44.beta2,
                                mpcar.eigen39.beta2))

data.mbym2 <- data.frame(prior=rep("MBYM2", 2100), model=model, 
                         beta=c(mbym2.beta2, mbym2.eigen64.beta2, mbym2.eigen59.beta2,
                                mbym2.eigen54.beta2, mbym2.eigen49.beta2, mbym2.eigen44.beta2,
                                mbym2.eigen39.beta2))



data.Scenario5 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario5$model <- factor(data.Scenario5$model, levels=c("M-Spatial", "M-SpatPlus64", 
                                                              "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", 
                                                              "M-SpatPlus39"))
data.Scenario5$prior <- factor(data.Scenario5$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario5$line <- c(rep(0.2, nrow(data.Scenario5)))
data.Scenario5$scenario <- rep("Scenario 5", nrow(data.Scenario5))




#########################################
# Merge datasets of different scenarios #
#########################################
data <- rbind(data.Scenario3, data.Scenario4, data.Scenario5)
data$model <- factor(data$model, levels=c("M-Spatial", "M-SpatPlus64", "M-SpatPlus59", 
                    "M-SpatPlus54", "M-SpatPlus49", "M-SpatPlus44", "M-SpatPlus39"))
data$prior <- factor(data$prior, levels=c("MICAR", "MPCAR", "MBYM2"))



######################################################
# Figure B5: boxplot of the estimated means of beta 1 #
######################################################
pdf(file="FigureB6.pdf", width=14, height=8)

p <- ggplot(data, aes(x = model, y = beta, fill=prior)) + 
  geom_boxplot(position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9),
        text = element_text(size = 20),
        panel.spacing = unit(1.3, "lines")) + 
  xlab("") + ylab(TeX('$\\beta_{2}$'))  +
  scale_y_continuous(limits=c(0.07, 0.48), breaks=c(0.1,0.2,0.3,0.4))

p + geom_hline(aes(yintercept = c(rep(0.2,nrow(data)))), color="#e31a1c", linetype="dashed") + 
  facet_wrap(~scenario, nrow=4, strip.position="right") + 
  theme(strip.background =element_rect(fill="lightskyblue")) + 
  scale_fill_manual(values = c("#fbb4ae","#b3cde3","#ccebc5")) 


dev.off()

