
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load packages
library(INLA)
library(ggplot2)
library(latex2exp)


n.sim <- 300 


beta2 <- function(x){
  data.frame(alpha.mean=x$summary.fixed[4,1])
}



##############
# Scenario 1 #
##############

# MICAR
load("results_Scenario1/MICAR_SexRatio_2df.Rdata")
load("results_Scenario1/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta2 <- do.call(rbind, lapply(micar.bart, beta2))[, 1]
micar.eigen64.beta2 <- do.call(rbind, lapply(micar.bart.eigen64, beta2))[, 1]
micar.eigen59.beta2 <- do.call(rbind, lapply(micar.bart.eigen59, beta2))[, 1]
micar.eigen54.beta2 <- do.call(rbind, lapply(micar.bart.eigen54, beta2))[, 1]
micar.eigen49.beta2 <- do.call(rbind, lapply(micar.bart.eigen49, beta2))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, micar.bart.eigen49)


# MPCAR
load("results_Scenario1/MPCAR_SexRatio_2df.Rdata")
load("results_Scenario1/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta2 <- do.call(rbind, lapply(mpcar.bart, beta2))[, 1]
mpcar.eigen64.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta2))[, 1]
mpcar.eigen59.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta2))[, 1]
mpcar.eigen54.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta2))[, 1]
mpcar.eigen49.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta2))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, mpcar.bart.eigen49)


# MBYM2
load("results_Scenario1/MBYM2_SexRatio_2df.Rdata")
load("results_Scenario1/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta2 <- do.call(rbind, lapply(mbym2.bart, beta2))[, 1]
mbym2.eigen64.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta2))[, 1]
mbym2.eigen59.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta2))[, 1]
mbym2.eigen54.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta2))[, 1]
mbym2.eigen49.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta2))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, mbym2.bart.eigen49)



###################################
# Scenario 1: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 1500), model=model, 
                         beta=c(micar.beta2, micar.eigen64.beta2, micar.eigen59.beta2,
                                micar.eigen54.beta2, micar.eigen49.beta2))

data.mpcar <- data.frame(prior=rep("MPCAR", 1500), model=model, 
                         beta=c(mpcar.beta2, mpcar.eigen64.beta2, mpcar.eigen59.beta2,
                                mpcar.eigen54.beta2, mpcar.eigen49.beta2))

data.mbym2 <- data.frame(prior=rep("MBYM2", 1500), model=model, 
                         beta=c(mbym2.beta2, mbym2.eigen64.beta2, mbym2.eigen59.beta2,
                                mbym2.eigen54.beta2, mbym2.eigen49.beta2))



data.Scenario1 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario1$model <- factor(data.Scenario1$model, levels=c("M-Spatial", "M-SpatPlus64", 
                                "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49"))
data.Scenario1$prior <- factor(data.Scenario1$prior, levels=c("MICAR", "MPCAR", "MBYM2"))

data.Scenario1$line <- c(rep(-0.2,nrow(data.Scenario1)))
data.Scenario1$scenario <- rep("Scenario 1", nrow(data.Scenario1))



##############
# Scenario 2 #
##############

# MICAR
load("results_Scenario2/MICAR_SexRatio_2df.Rdata")
load("results_Scenario2/MICAR_SpatPlus_SexRatio_2df.Rdata")

micar.beta2 <- do.call(rbind, lapply(micar.bart, beta2))[, 1]
micar.eigen64.beta2 <- do.call(rbind, lapply(micar.bart.eigen64, beta2))[, 1]
micar.eigen59.beta2 <- do.call(rbind, lapply(micar.bart.eigen59, beta2))[, 1]
micar.eigen54.beta2 <- do.call(rbind, lapply(micar.bart.eigen54, beta2))[, 1]
micar.eigen49.beta2 <- do.call(rbind, lapply(micar.bart.eigen49, beta2))[, 1]

rm(micar.bart, micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, micar.bart.eigen49)


# MPCAR
load("results_Scenario2/MPCAR_SexRatio_2df.Rdata")
load("results_Scenario2/MPCAR_SpatPlus_SexRatio_2df.Rdata")

mpcar.beta2 <- do.call(rbind, lapply(mpcar.bart, beta2))[, 1]
mpcar.eigen64.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen64, beta2))[, 1]
mpcar.eigen59.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen59, beta2))[, 1]
mpcar.eigen54.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen54, beta2))[, 1]
mpcar.eigen49.beta2 <- do.call(rbind, lapply(mpcar.bart.eigen49, beta2))[, 1]

rm(mpcar.bart, mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, mpcar.bart.eigen49)


# MBYM2
load("results_Scenario2/MBYM2_SexRatio_2df.Rdata")
load("results_Scenario2/MBYM2_SpatPlus_SexRatio_2df.Rdata")

mbym2.beta2 <- do.call(rbind, lapply(mbym2.bart, beta2))[, 1]
mbym2.eigen64.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen64, beta2))[, 1]
mbym2.eigen59.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen59, beta2))[, 1]
mbym2.eigen54.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen54, beta2))[, 1]
mbym2.eigen49.beta2 <- do.call(rbind, lapply(mbym2.bart.eigen49, beta2))[, 1]

rm(mbym2.bart, mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, mbym2.bart.eigen49)



###################################
# Scenario 2: prepare a dataframe #
###################################
model <- c(rep("M-Spatial", n.sim), rep("M-SpatPlus64", n.sim), rep("M-SpatPlus59", n.sim),
           rep("M-SpatPlus54", n.sim), rep("M-SpatPlus49", n.sim))


data.micar <- data.frame(prior=rep("MICAR", 1500), model=model, 
                         beta=c(micar.beta2, micar.eigen64.beta2, micar.eigen59.beta2,
                                micar.eigen54.beta2, micar.eigen49.beta2))

data.mpcar <- data.frame(prior=rep("MPCAR", 1500), model=model, 
                         beta=c(mpcar.beta2, mpcar.eigen64.beta2, mpcar.eigen59.beta2,
                                mpcar.eigen54.beta2, mpcar.eigen49.beta2))

data.mbym2 <- data.frame(prior=rep("MBYM2", 1500), model=model, 
                         beta=c(mbym2.beta2, mbym2.eigen64.beta2, mbym2.eigen59.beta2,
                                mbym2.eigen54.beta2, mbym2.eigen49.beta2))



data.Scenario2 <- rbind(data.micar, data.mpcar, data.mbym2)
data.Scenario2$model <- factor(data.Scenario2$model, levels=c("M-Spatial", "M-SpatPlus64", 
                               "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49"))
data.Scenario2$prior <- factor(data.Scenario2$prior, levels=c("MICAR", "MPCAR", "MBYM2"))
data.Scenario2$line <- c(rep(-0.2,nrow(data.Scenario2)))
data.Scenario2$scenario <- rep("Scenario 2", nrow(data.Scenario2))




#########################################
# Merge datasets of different scenarios #
#########################################
data <- rbind(data.Scenario1, data.Scenario2)
data$model <- factor(data$model, levels=c("M-Spatial", "M-SpatPlus64", 
                     "M-SpatPlus59", "M-SpatPlus54", "M-SpatPlus49"))
data$prior <- factor(data$prior, levels=c("MICAR", "MPCAR", "MBYM2"))



######################################################
# Figure 5: boxplot of the estimated means of beta 1 #
######################################################
pdf(file="Figure5.pdf", width=14, height=8)

p <- ggplot(data, aes(x = model, y = beta, fill=prior)) + 
  geom_boxplot(position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9),
        text = element_text(size = 20),
        panel.spacing = unit(1.3, "lines")) + 
  xlab("") + ylab(TeX('$\\beta_{2}$'))  +
  scale_y_continuous(limits=c(-0.5, -0.04), breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1))

p + geom_hline(aes(yintercept = c(rep(-0.2,nrow(data)))), color="#e31a1c", linetype="dashed") + 
  facet_wrap(~scenario, nrow=2, strip.position="right") + 
  theme(strip.background =element_rect(fill="lightskyblue")) + 
  scale_fill_manual(values = c("#fbb4ae","#b3cde3","#ccebc5")) 

dev.off()

