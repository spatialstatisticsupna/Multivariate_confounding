
###########################################
# Mmodels - BYM2 (BARTLETT DECOMPOSITION) #
###########################################
Mmodel.bym2.bartlett <- function(
  cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"), 
      theta=NULL){
  
  envir <- parent.env(environment())
  if(!exists("cache.done", envir=envir)){
    Qi <- Diagonal(x=colSums(W))-W
    Ue <- eigen(Qi)$vectors
    eval <- exp(mean(log(diag(INLA:::inla.ginv(Qi)))))
    eigen <- eigen(Qi)$values[-nrow(Qi)] 
    
    assign("Ue", Ue, envir = envir)
    assign("eval", eval, envir = envir)
    assign("eigen", eigen, envir = envir)
    assign("cache.done", TRUE, envir = envir)
  }

  # ------------------------------------------------------------- #
  # Covert the hyperparameters from internal scale to model scale #
  # ------------------------------------------------------------- #
  interpret.theta <- function(){
    alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-theta[as.integer(1:J)]))
     
    diag.N <- sapply(theta[J+1:J], function(x) {exp(x)})
    no.diag.N <- theta[2*J+1:(J*(J-1)/2)]
    
    N <- diag(diag.N) 
    N[lower.tri(N)] <- no.diag.N
    
    Covar <- N %*% t(N)
    
    e <- eigen(Covar)
    M <- t(e$vectors %*% diag(sqrt(e$values)))

    return(list(alpha=alpha, Covar=Covar, M=M))
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
    M.inv <- solve(param$M)
    MI <- kronecker(M.inv, Ue)
    BlockIW <- 
      Matrix::bdiag(lapply(1:J, function(i) {
        Matrix::Diagonal(x=c(eigen/(eigen+param$alpha[i]*(eval-eigen)), 1/(1-param$alpha[i])))
      }))
    Q <- (MI %*% BlockIW) %*% t(MI)
    Q <- inla.as.sparse(Q)
    return (Q)
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
    # Uniform prior in (alpha.min, alpha.max) on model scale
    val <-  sum(-theta[1:J] - 2*log(1+exp(-theta[1:J])))
    # n^2_jj ~ chisq(J-j+1) (J degrees of freedom)
    val <- val + J*log(2) + 2*sum(theta[J+1:J]) + sum(dchisq(exp(2*theta[J+1:J]), df=J-1:J+1, log=TRUE))
    # n_ji ~ N(0,1)
    val <- val + sum(dnorm(theta[(2*J)+1:(J*(J-1)/2)], mean=0, sd=1, log=TRUE))
    return(val)
  }

  # --------------------- #
  # Return initial values #
  # --------------------- #
  initial <- function(){
    p <- (0.9-alpha.min)/(alpha.max-alpha.min)
    
    return(c(rep(log(p/(1-p)),J), as.vector(initial.values)))
  }
  
  # --------------- #
  # Function quit() #
  # --------------- #
  quit <- function(){ return(invisible()) }
  
  if(!length(theta)) theta <- initial()
  val <- do.call(match.arg(cmd), args=list())
  
  return(val)
}



