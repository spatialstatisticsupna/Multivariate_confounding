
rm(list=ls())
setwd("")

# Load packages
library(sf)
library(RColorBrewer)
library(tmap)


# Define the scenario 
Scenario <- 3



#################
# Load the data #
#################
load(paste0("Simulated_data/SimuStudy1_Scenario", Scenario, ".Rdata"))
head(data)

S <- length(unique(data$dist))



###############################################
# Figure 2: plot the covariates X1, X2 and X3 #
###############################################
carto_UP$X1 <- data$X1[1:70]
carto_UP$X2 <- data$X[1:S]
carto_UP$X3 <- data$X[(S+1):(2*S)]


paleta <- brewer.pal(9,"YlOrRd")
values <- c(-2.37,-0.8,-0.5,-0.25,0,0.25,0.5,0.8, 1.5,2.66)


tmap_mode("plot")
plot <- tm_shape(carto_UP) +
  tm_polygons(col=c("X1", "X2", "X3"), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left", border.col="black") +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1,
            panel.labels=c("X1","X2", "X3"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F, outer.margins=0) +
  tm_facets(ncol=3, nrow=1)
print(plot)

tmap_save(plot, file="FigureA1a.pdf", width=12, height=4)
# tmap_save(plot, file="FigureA1b.pdf", width=12, height=4)
# tmap_save(plot, file="FigureA1c.pdf", width=12, height=4)
# tmap_save(plot, file="FigureA1d.pdf", width=12, height=4)


