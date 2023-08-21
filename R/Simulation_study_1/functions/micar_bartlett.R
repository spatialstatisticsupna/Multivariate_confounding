
###########################################
# Mmodels - ICAR (BARTLETT DECOMPOSITION) #
###########################################
Mmodel.icar.bartlett <- function(
  cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"), 
  theta=NULL){
  
  envir <- parent.env(environment())
  if(!exists("cache.done", envir=envir)){
    IW <- Diagonal(x=colSums(W))-W
    assign("IW", IW, envir=envir)
    assign("cache.done", TRUE, envir=envir)
  }
  
  # ------------------------------------------------------------- #
  # Covert the hyperparameters from internal scale to model scale #
  # ------------------------------------------------------------- #
  interpret.theta <- function(){
    diag.N <- sapply(theta[1:J], function(x) {exp(x)})
    no.diag.N <- theta[(J+1):((J*(J+1))/2)]
    
    N <- diag(diag.N) 
    N[lower.tri(N)] <- no.diag.N
    
    Covar <- N %*% t(N)
    
    e <- eigen(Covar)
    M <- t(e$vectors %*% diag(sqrt(e$values)))

    return(list(Covar=Covar, M=M))
  }
  
  # ---------------------------------------------------------------------------- #
  # Graph of the precision matrix, i.e. a 0/1 representation of precision matrix #
  # ---------------------------------------------------------------------------- #
  graph <- function(){ return(Q()) }
  
  # ---------------- #
  # Precision matrix #
  # ---------------- #
  Q <- function(){
   param <- interpret.theta()
   Covar.inv <- solve(param$Covar)
   Q <- kronecker(Covar.inv, IW)
   return(Q)
  }
  
  # ------------- #
  # Mean of model #
  # ------------- #
  mu <- function(){ return(numeric(0)) }
  
  # -------------- #
  # log.norm.const #
  # -------------- #
  log.norm.const <- function(){
    val <- numeric(0)
    return(val)
  }
  
  # -------------------------------------------------------------- #
  # Define the log-prior for the hyperparameters in internal scale #
  # -------------------------------------------------------------- #
  log.prior <- function(){
    param <- interpret.theta()
    # n^2_jj ~ chisq(J-j+1) (J degrees of freedom)
    val <- sum(dchisq(exp(2*theta[1:J]), df=as.integer(J-1:J+1), log=TRUE)) + J*log(2) + 2*sum(theta[1:J])
    # n_jk ~ N(0,1)
    val <- val + sum(dnorm(theta[(J+1):((J*(J+1))/2)], mean=0, sd=1, log=TRUE))
    return(val)
  }

  # --------------------- #
  # Return initial values #
  # --------------------- #
  initial <- function(){ return(as.vector(initial.values)) }
  
  # --------------- #
  # Function quit() #
  # --------------- #
  quit <- function(){ return(invisible()) }
  
  if(!length(theta)) theta <- initial()
  val <- do.call(match.arg(cmd), args=list())
  
  return(val)
}


