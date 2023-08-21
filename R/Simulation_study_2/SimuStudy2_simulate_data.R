
################################################################################
#                    Simulation study 2: data simulation                       #
################################################################################

rm(list=ls())
setwd("")

# Load packages
library(MASS)
library(spdep)



#################
# Load the data #
#################
load("Data/data_UttarPradesh_2011.RData")

J <- length(unique(data$Crime))
S <- length(unique(data$dist))



###################################
# Define spatial structure matrix #
###################################
nb.mat <- nb2mat(poly2nb(carto_UP), style ="B")
nbs <- unlist(lapply(poly2nb(carto_UP), function(x) length(x)))

Q.xi <- diag(nbs)-nb.mat



########################################################
# Define the variance-covariance matrix between crimes #
########################################################

# Choose the scenario you are interested in

# Scenario 1
cor.crimes <- 0.7

# Scenario 2, Scenario 4, Scenario 5
# cor.crimes <- 0.3

# Scenario 3
# cor.crimes <- 0.5


sd.1 <- sqrt(0.9)
sd.2 <- sqrt(0.2)

cor.mat <- matrix(c(1,cor.crimes,cor.crimes,1),2,2)


# Between crimes covariance matrix
Sigma.b <- diag(c(sd.1, sd.2))%*%cor.mat%*%diag(c(sd.1, sd.2))
# Within crimes covariance matrix
Sigma.w <- ginv(Q.xi)


Sigma <- kronecker(Sigma.b, Sigma.w)



################################
# Simulate the spatial effects #
################################
set.seed(1302)
xi <- mvrnorm(1, rep(0,2*S), Sigma)


cor(xi[1:S], xi[(S+1):(2*S)])



#############
# Functions #
#############
complement2 <- function(y, rho, x, threshold=1e-12) {
  if(!is.matrix(y)) y <- matrix(y, ncol=1)
  d <- ncol(y)
  n <- nrow(y)
  y <- scale(y, center=FALSE) # Makes computations simpler
  if (missing(x)) x <- rnorm(n)
  
  e <- residuals(lm(x ~ y))
  y.dual <- with(svd(y), (n-1)*u %*% diag(ifelse(d > threshold, 1/d, 0)) %*% t(v))
  sigma2 <- c((1 - rho %*% cov(y.dual) %*% rho) / var(e))
  
  if (sigma2 >= 0) {
    sigma <- sqrt(sigma2) 
    z <- y.dual %*% rho + sigma*e
  } else {
    warning("Correlations are impossible.")
    z <- rep(0, n)
  }
  return(z)
}



####################################################################
# Simulate X1* with a desired correlation with the spatial effects #
####################################################################

# Choose the scenario you are interested in

# Scenario 1, Scenario 3
X1 <- complement2(cbind(xi[1:S], xi[(S+1):(2*S)]), rho=c(0.5, 0.7))

# Scenario 2
# X1 <- complement2(cbind(xi[1:S], xi[(S+1):(2*S)]), rho=c(0.3, 0.5))

# Scenario 4
# X1 <- complement2(cbind(xi[1:S], xi[(S+1):(2*S)]), rho=c(0.3, 0.7))

# Scenario 5
# X1 <- complement2(cbind(xi[1:S], xi[(S+1):(2*S)]), rho=c(0, 0.7))


cor(X1, xi[1:S])
cor(X1, xi[(S+1):(2*S)])


# Add X1* and spatial effects to the dataset
data$xi <- xi
data$X1.aux <- c(X1, X1)



#########################
# Compute the log risks #
#########################
log.risk.crime1 <- 0.12 +0.15*data$X1.aux[1:S] + data$xi[1:S]
lambda.crime1 <-  data$exp[1:S]*exp(log.risk.crime1)


log.risk.crime2 <- 0.03 +0.2*data$X1.aux[1:S] + data$xi[(S+1):(2*S)]
lambda.crime2 <-  data$exp[(S+1):(2*S)]*exp(log.risk.crime2)


lambda <- c(lambda.crime1, lambda.crime2)
log.risk <- c(log.risk.crime1, log.risk.crime2)



#######################
# Simulate the counts # 
#######################
n.sim <- 300
simu.O <- vector("list", n.sim)
simu.O.crime1 <- vector("list", n.sim)
simu.O.crime2 <- vector("list", n.sim)

for(i in 1:n.sim) {
  set.seed(10+i)
  O <- rpois(2*S, lambda)
  simu.O[[i]] <- O
  simu.O.crime1[[i]] <- O[1:S]
  simu.O.crime2[[i]] <- O[(S+1):(2*S)]
}


# Folder to save results
if(!file.exists("Simulated_data")) dir.create("Simulated_data")


# Define the scenario
Scenario <- 1


# Save simulated data
save(simu.O, simu.O.crime1, simu.O.crime2, log.risk.crime1, log.risk.crime2, log.risk, data, carto_UP, 
     file=paste0("Simulated_data/SimuStudy2_Scenario", Scenario, ".Rdata"))



