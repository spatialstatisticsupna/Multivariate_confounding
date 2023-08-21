
################################################################################
#                    Simulation study 1: data simulation                       #
################################################################################

rm(list=ls())
setwd("")

# Load packages
library(MASS)



#################
# Load the data #
#################
load("Data/data_UttarPradesh_2011.RData")

J <- length(unique(data$Crime))
S <- length(unique(data$dist))


X1 <- data$X1[1:S]



#############
# Functions #
#############
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

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

orth <- function(X,n.vectors=1){
  n <- nrow(X)
  P <- diag(n)-X%*%solve(t(X)%*%X)%*%t(X)
  
  res <- P%*%matrix(rnorm(n*n.vectors),n,n.vectors)
  return(res)
}



######################
# Simulate X2 and X3 #
######################

# Choose the scenario you are interested in

# Scenario 1: define the correlations
cor.X1.X2 <- 0.5
cor.X1.X3 <- 0.7
cor.X2.X3 <- 0.7

# Scenario 2: define the correlations
# cor.X1.X2 <- 0.3
# cor.X1.X3 <- 0.5
# cor.X2.X3 <- 0.3

# Scenario 3: define the correlations
# cor.X1.X2 <- 0.5
# cor.X1.X3 <- 0.7
# cor.X2.X3 <- 0.5

# Scenario 4: define the correlations
# cor.X1.X2 <- 0.3
# cor.X1.X3 <- 0.7
# cor.X2.X3 <- 0.3

# Scenario 5: define the correlations
# cor.X1.X2 <- 0
# cor.X1.X3 <- 0.7
# cor.X2.X3 <- 0.7

# Scenario 6: define the correlations
# cor.X1.X2 <- 0
# cor.X1.X3 <- 0.7
# cor.X2.X3 <- 0.3



set.seed(1289)

# Generate X3
X3 <- scale(complement(X1,rho=cor.X1.X3))
X <- cbind(X1,X3)
cor(X)

# Generate X2
X2 <- complement2(cbind(X1,X3),rho=c(cor.X1.X2, cor.X2.X3))
X <- cbind(X1,X2,X3)
cor(X)


# Add X2 and X3 to the dataset
data$X <- c(X2, X3)



#########################
# Compute the log risks #
#########################
log.risk.crime1 <- -0.12 - 0.15*X1 - 0.3*X2
lambda.crime1 <-  data$exp[1:S]*exp(log.risk.crime1)


log.risk.crime2 <- -0.03 - 0.2*X1 - 0.3*X3
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
     file=paste0("Simulated_data/SimuStudy1_Scenario", Scenario, ".Rdata"))




