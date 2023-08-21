
rm(list=ls())
setwd("")

n.sim <- 300

# Define the scenario
Scenario <- 1



###################
# MICAR M-Spatial #
###################
micar.bart <- vector("list", n.sim)

for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MICAR_Spat_", i, ".Rdata"))
  micar.bart[[i]] <- micar.spat
}

save(micar.bart, file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MICAR_SexRatio_2df.Rdata"))



####################
# MICAR M-SpatPlus #
####################
micar.bart.eigen64 <- vector("list", n.sim)
micar.bart.eigen59 <- vector("list", n.sim)
micar.bart.eigen54 <- vector("list", n.sim)
micar.bart.eigen49 <- vector("list", n.sim)
micar.bart.eigen44 <- vector("list", n.sim)
micar.bart.eigen39 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MICAR_SpatPlus_", i, ".Rdata"))
  micar.bart.eigen64[[i]] <- micar.eigen64
  micar.bart.eigen59[[i]] <- micar.eigen59
  micar.bart.eigen54[[i]] <- micar.eigen54
  micar.bart.eigen49[[i]] <- micar.eigen49
  micar.bart.eigen44[[i]] <- micar.eigen44
  micar.bart.eigen39[[i]] <- micar.eigen39
}

save(micar.bart.eigen64, micar.bart.eigen59, micar.bart.eigen54, micar.bart.eigen49,
     micar.bart.eigen44, micar.bart.eigen39, 
     file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MICAR_SpatPlus_SexRatio_2df.Rdata"))



###################
# MPCAR M-Spatial #
###################
mpcar.bart <- vector("list", n.sim)

for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MPCAR_Spat_", i, ".Rdata"))
  mpcar.bart[[i]] <- mpcar.spat
}

save(mpcar.bart, file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MPCAR_SexRatio_2df.Rdata"))



####################
# MPCAR M-SpatPlus #
####################
mpcar.bart.eigen64 <- vector("list", n.sim)
mpcar.bart.eigen59 <- vector("list", n.sim)
mpcar.bart.eigen54 <- vector("list", n.sim)
mpcar.bart.eigen49 <- vector("list", n.sim)
mpcar.bart.eigen44 <- vector("list", n.sim)
mpcar.bart.eigen39 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MPCAR_SpatPlus_", i, ".Rdata"))
  mpcar.bart.eigen64[[i]] <- mpcar.eigen64
  mpcar.bart.eigen59[[i]] <- mpcar.eigen59
  mpcar.bart.eigen54[[i]] <- mpcar.eigen54
  mpcar.bart.eigen49[[i]] <- mpcar.eigen49
  mpcar.bart.eigen44[[i]] <- mpcar.eigen44
  mpcar.bart.eigen39[[i]] <- mpcar.eigen39
}


save(mpcar.bart.eigen64, mpcar.bart.eigen59, mpcar.bart.eigen54, mpcar.bart.eigen49, 
     mpcar.bart.eigen44, mpcar.bart.eigen39,
     file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MPCAR_SpatPlus_SexRatio_2df.Rdata"))



###################
# MBYM2 M-Spatial #
###################
mbym2.bart <- vector("list", n.sim)

for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MBYM2_Spat_", i, ".Rdata"))
  mbym2.bart[[i]] <- mbym2.spat
}

save(mbym2.bart, file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MBYM2_SexRatio_2df.Rdata"))



####################
# MBYM2 M-SpatPlus #
####################
mbym2.bart.eigen64 <- vector("list", n.sim)
mbym2.bart.eigen59 <- vector("list", n.sim)
mbym2.bart.eigen54 <- vector("list", n.sim)
mbym2.bart.eigen49 <- vector("list", n.sim)
mbym2.bart.eigen44 <- vector("list", n.sim)
mbym2.bart.eigen39 <- vector("list", n.sim)


for (i in 1:n.sim){
  print(i)
  load(paste0("Simulation_study_2/results_Scenario", Scenario, "/MBYM2_SpatPlus_", i, ".Rdata"))
  mbym2.bart.eigen64[[i]] <- mbym2.eigen64
  mbym2.bart.eigen59[[i]] <- mbym2.eigen59
  mbym2.bart.eigen54[[i]] <- mbym2.eigen54
  mbym2.bart.eigen49[[i]] <- mbym2.eigen49
  mbym2.bart.eigen44[[i]] <- mbym2.eigen44
  mbym2.bart.eigen39[[i]] <- mbym2.eigen39
}

save(mbym2.bart.eigen64, mbym2.bart.eigen59, mbym2.bart.eigen54, mbym2.bart.eigen49, 
     mbym2.bart.eigen44, mbym2.bart.eigen39,
     file=paste0("Simulation_study_2/results_Scenario", Scenario, "/MBYM2_SpatPlus_SexRatio_2df.Rdata"))






