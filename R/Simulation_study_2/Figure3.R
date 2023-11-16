
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(sf)
library(RColorBrewer)
library(tmap)
library(latex2exp)


# Define the scenario
scenario <- "Scenario1"
# scenario <- "Scenario2"



#################
# Load the data #
#################
load(paste0("Simulated_data/SimuStudy2_", scenario, ".Rdata"))
head(data)

S <- length(unique(data$dist))



##################################################################
# Figure 3: plot simulated covariate X1* and the spatial effects #
##################################################################
carto_UP$X1.aux <- data$X1.aux[1:S]
carto_UP$xi.crime1 <- data$xi[1:S]
carto_UP$xi.crime2 <- data$xi[(S+1):(2*S)]



paleta <- brewer.pal(9,"YlOrRd")
values <- c(-3.21,-0.8,-0.5,-0.25,0,0.25,0.5,0.8, 1.5,3.27)


tmap_mode("plot")
plot <- tm_shape(carto_UP) +
  tm_polygons(col=c("X1.aux", "xi.crime1", "xi.crime2"), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left", border.col="black") +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1,
            panel.labels=c("X1*",TeX('$\\bf{Theta}_{1}$'), TeX('$\\Theta_{2}$')), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F, outer.margins=0) +
  tm_facets(ncol=3, nrow=1)

print(plot)

if (scenario=="Scenario1"){
  tmap_save(plot, file="Figure3a.pdf", width=12, height=4)
}
if (scenario=="Scenario2"){
  tmap_save(plot, file="Figure3b.pdf", width=12, height=4)
}

