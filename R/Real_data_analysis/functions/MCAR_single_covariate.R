
options(dplyr.summarise.inform = FALSE)


##########################################################################
# Fit a spatial multivariate Poisson mixed model to areal count data     #
# where dependence between spatial patterns of the diseases is addressed #
# through the use of M-models (Botella-Rocamora et al, 2015)             #
##########################################################################

MCAR.Covariate <- function(carto=NULL, data=NULL, ID.area=NULL, ID.disease=NULL,
                             O=NULL, E=NULL, spatial="intrinsic.bartlett", 
                             strategy="simplified.laplace", covariate=NULL){
  
  # Check for errors #
  if(is.null(carto))
    stop("the 'carto' argument is missing")
  if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
    stop("the 'carto' argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
  if(is.null(ID.area))
    stop("the 'ID.area' argument is missing")
  if(is.null(ID.disease))
    stop("the 'ID.disease' argument is missing")
  if(is.null(O))
    stop("the 'O' argument is missing")
  if(is.null(E))
    stop("the 'E' argument is missing")
  if(!(spatial %in% c("intrinsic.bartlett", "proper.bartlett","bym2.bartlett")))
    stop("invalid 'spatial' argument")
  if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
    stop("invalid 'strategy' argument")
  
  if(!ID.area %in% colnames(carto))
    stop(sprintf("'%s' variable not found in carto object", ID.area))
  if(!ID.area %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object", ID.area))
  if(!ID.disease %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object", ID.disease))
  if(!O %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object", O))
  if(!E %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object", E))
  
  
  cat("STEP 1: Pre-processing data\n")
  
  # Transform 'SpatialPolygonsDataFrame' object to 'sf' class #
  carto <- st_as_sf(carto)
  
  # Order the data #
  data[, ID.disease] <- paste(sprintf("%02d", as.numeric(as.character(data[, ID.disease]))))
  data <- data[order(data[, ID.disease], data[, ID.area]), ]
  
  # Define spatial adjacency matrix #
  carto.nb <- poly2nb(carto)
  Ws <- as(nb2listw(carto.nb, style="B"),"CsparseMatrix")
  
  # Define S and J #
  S <- length(carto.nb)
  J <- length(unique(data[,ID.disease]))

  # Define hyperprior distributions #
  sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"
  
  lunif = "expression:
          a = 1;
          b = 1;
          beta = exp(theta)/(1+exp(theta));
          logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
          log_jacobian = log(beta*(1-beta));
          return(logdens+log_jacobian)"

  # Formula for INLA model #
  form <- "O ~ -1+"
  form <- paste(form, paste(paste0("I", 1:J), collapse="+"), sep="")
  form <- paste(form, "+ ")
  form <- paste(form, paste(paste0("X.crime", 1:J), collapse="+"), sep="")
  form <- paste(form, "+ f(idx, model=Mmodel.s, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,J)))")
  formula <- stats::as.formula(form)
  
  
  cat("STEP 2: Fitting the model with INLA (this may take a while...)\n")

  # Identifiability constraints for each disease #
  A.constr.s <- kronecker(diag(J), matrix(1,1,S))
  
  # Define data frame for INLA model #
  data.INLA <- data.frame(O=data[,O], E=data[,E], Area=data[,ID.area], Disease=data[,ID.disease],
                          ID.area=rep(1:S, J), ID.disease=rep(1:J, each=S), X=data[, covariate])
  
  intercepts <- dummy_cols(data.INLA$ID.disease)[,-1]
  intercepts[intercepts==0] <- NA
  colnames(intercepts) <- paste0("I", 1:J)
  data.INLA <- cbind(data.INLA, intercepts)
  
  data.INLA$idx <- (data.INLA$ID.disease-1)*S + data.INLA$ID.area
  
  X <- intercepts*data.INLA$X
  colnames(X) <- paste0("X.crime", 1:J)
  data.INLA <- cbind(data.INLA, X)

  # Initial values for spatial M-model #
  SMR <- data[,O]/data[,E]
  Sigma <- cov(matrix(SMR,S,J,byrow=F))
  N <- t(chol(Sigma))
  Rho <- cor(matrix(SMR,S,J,byrow=F))
  
  # Define selected M-model #
  if(spatial=="intrinsic.bartlett"){
    initial.values.s <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))
    Mmodel.s <- inla.rgeneric.define(Mmodel.icar.bartlett, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s)
  }
  if(spatial=="proper.bartlett"){
    initial.values.s <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))
    Mmodel.s <- inla.rgeneric.define(Mmodel.pcar.bartlett, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
  }
  if(spatial=="bym2.bartlett"){
    initial.values.s <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))
    Mmodel.s <- inla.rgeneric.define(Mmodel.bym2.bartlett, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
  }
  
  # Fit the INLA model  #
  Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                control.inla=list(strategy=strategy))
  
  Model$Mmodel <- list(spatial=spatial)
  
  # Compute spatial between-disease correlations and marginal variances #
  if(spatial %in% c("intrinsic.bartlett", "proper.bartlett", "bym2.bartlett")){
    Mmodel.compute <- compute.cor(Model, n.sample=1000)
    Model$summary.cor <- Mmodel.compute$summary.cor
    Model$marginals.cor <- Mmodel.compute$marginals.cor
    Model$summary.var <- Mmodel.compute$summary.var
    Model$marginals.var <- Mmodel.compute$marginals.var
  }
  
  return(Model)
}


#######################
# Auxiliary functions #
#######################

tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


compute.summary <- function(marginal){
  aux <- data.frame(inla.emarginal(function(x) x, marginal),
                    sqrt(inla.emarginal(function(x) x^2, marginal)-inla.emarginal(function(x) x, marginal)^2),
                    inla.qmarginal(0.025,marginal),
                    inla.qmarginal(0.5,marginal),
                    inla.qmarginal(0.975,marginal),
                    inla.mmarginal(marginal),
                    1-inla.pmarginal(0,marginal))
  colnames(aux) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode","1 cdf")
  aux
}


compute.cor <- function(model, n.sample=10000){
  o <- tryCatch.W.E({
    J <- length(unique(model$.args$data$ID.disease))
    hyperpar.sample <- inla.hyperpar.sample(n.sample, model, improve.marginals=TRUE)
    hyperpar.sample.s <- hyperpar.sample[, grep("idx", colnames(hyperpar.sample))]
    
    # Covariance/correlation matrix for spatial M-model #
    if(model$Mmodel$spatial=="intrinsic.bartlett"){
      hyperpar.sample.s[,1:J] <- exp(hyperpar.sample.s[,1:J])
      hyperpar.sample.s <- split(hyperpar.sample.s[,seq(J*(J+1)/2)], seq(nrow(hyperpar.sample.s)))
    }
    if(model$Mmodel$spatial %in% c("proper.bartlett", "bym2.bartlett")){
      hyperpar.sample.s[,seq(J+1,2*J)]<- exp(hyperpar.sample.s[,seq(J+1,2*J)])
      hyperpar.sample.s <- split(hyperpar.sample.s[,seq(1+J,J+J*(J+1)/2)], seq(nrow(hyperpar.sample.s)))
    }

    param.sample.s <- lapply(hyperpar.sample.s, function(x){
      N <- diag(x[seq(J)]) 
      N[lower.tri(N, diag=FALSE)] <- x[-seq(J)]
      Sigma <- N %*% t(N)
      Rho <- cov2cor(Sigma)
      Rho.values <- Rho[lower.tri(Rho)]
      
      return(list(sigma=diag(Sigma),rho=Rho.values))
    })
    
    cor.sample.s <- do.call(rbind,lapply(param.sample.s, function(x) x$rho))
    cor.density.s <- apply(cor.sample.s, 2, function(x) density(x, n=75, bw="SJ", from=-1, to=1))
      
    marginals.cor.s <- lapply(cor.density.s, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.cor.s) <- paste("rho.s",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep="")
    summary.cor.s <- do.call(rbind,lapply(marginals.cor.s, compute.summary))[,1:6]
      
    var.sample.s <- do.call(rbind,lapply(param.sample.s, function(x) x$sigma))
    var.density.s <- apply(var.sample.s, 2, function(x) density(x, n=75, bw="SJ", from=0))
      
    marginals.var.s <- lapply(var.density.s, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.var.s) <- paste("var.s",1:J,sep="")
    summary.var.s <- do.call(rbind,lapply(marginals.var.s, compute.summary))[,1:6]
  })
    
  Mmodel.compute <- list(summary.cor=summary.cor.s,
                         marginals.cor=marginals.cor.s,
                         summary.var=summary.var.s,
                         marginals.var=marginals.var.s)
  
  return(Mmodel.compute)
}
