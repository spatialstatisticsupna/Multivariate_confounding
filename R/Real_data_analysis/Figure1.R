
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(sf)
library(RColorBrewer)
library(tmap)



#################
# Load the data #
#################

# Crime: 1=rapes, 2=dowry deaths
load("../../Data/data_UttarPradesh_2011.Rdata")
head(data)

J <- length(unique(data$Crime))
S <- length(unique(data$dist))


carto_UP$X1 <- data$X1[1:S]



#############################################
# Compute the SIR of rapes and dowry deaths #
#############################################
carto_UP$SIR.rape <- data$obs[1:S]/data$exp[1:S]
carto_UP$SIR.dowry <- data$obs[(S+1):(2*S)]/data$exp[(S+1):(2*S)]



#######################################################################
# Figure 1: SIR of rapes, SIR of dowry deaths and Sex ratio covariate #
#######################################################################

# SIR of rapes and dowry deaths
summary(carto_UP[, c("SIR.rape", "SIR.dowry")])


paleta <- brewer.pal(9,"YlOrRd")
values1 <- c(0.12,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.37)


SIR.rape <- tm_shape(carto_UP) +
  tm_polygons(col=c("SIR.rape"), id="ID_area", 
              palette=paleta, border.alpha=1, title="", leyend.show=T, 
              legend.reverse=T, style="fixed", breaks=values1, interval.closure="left", 
              border.col="black") +
  tm_layout(main.title.position="left", legend.text.size=1,
            panel.labels=c("SIR rape"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F, 
            outer.margins=0) +
  tm_facets(ncol=1, nrow=1)


SIR.dowry <- tm_shape(carto_UP) +
  tm_polygons(col=c("SIR.dowry"), id="ID_area", 
              palette=paleta, border.alpha=1, title="",leyend.show=T, 
              legend.reverse=T, style="fixed", breaks=values1, interval.closure="left", 
              border.col="black") +
  tm_layout(main.title.position="left", legend.text.size=1,
            panel.labels=c("SIR dowry deaths"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F, 
            outer.margins=0) +
  tm_facets(ncol=1, nrow=1)



# Sex ratio
summary(carto_UP$X1)

values2 <- c(-1.36,-0.8,-0.5,-0.25,0,0.25,0.5,0.8, 1.5,2.66)


SexRatio <- tm_shape(carto_UP) +
  tm_polygons(col=c("X1"), 
              palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values2, interval.closure="left", 
              border.col="black") +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1,
            panel.labels=c("Sex ratio"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F, 
            outer.margins=0) +
  tm_facets(ncol=1, nrow=1)



# Merge all the figures
plot <- tmap_arrange(SIR.rape, SIR.dowry, SexRatio, ncol=3, nrow=1, outer.margins = 0) 
print(plot)


# Save the figure
tmap_save(plot, width=12, height=4, file="Figure1.pdf")





