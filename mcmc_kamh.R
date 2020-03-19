## Description
#   Implementation of the Kernel Adaptive Metropolis-Hastings (KAMH) algorithm proposed
#   by Sejdinovic et al. 2014.
#   [Reference: D. Sejdinovic, H. Strathmann, M.G. Lomeli, C. Andrieu, and A. Gretton. 
#   Kernel Adaptive Metropolis-Hastings. In International Conference on Machine Learning, 
#   JMLR W&CP 32(2), pages 1665--1673, 2014.]
## Arguments
#   (i)     target: function of the target pdf
#   (ii)    dimension: integer specifying the dimensionality of the target distribution
#   (iii)   X0: numerical vector specifying the initial state for the MCMC algorithm
#   (iv)    length.out: integer specifying output chain length, including the initial
#           state taken from the random walk M-H algorithm. Note that, this is not
#           necessarily the number of iterations, as this specifies chain length
#           returned after accounting for thinning.
#   (v)     max.subchain.length: integer specifying the maximum length of the sampled 
#           subchain used during the adaptation step
#   (vi)    gam: initial value for the scaling parameter gamma, part of the proposal's 
#           covariance matrix
#   (vii)   nu: scaling parameter nu, part of the proposal's covariance matrix
#   (viii)  nu.adapt: function which takes iteration number as an argument and
#           returns the value of a scaling parameter for adapting the nu parameter
#   (ix)    alpha: target acceptance rate used in adapting the nu parameter
#   (x)     adapt.prob: function which takes iteration number as an argument and 
#           returns numerical adaptation probability (requires that at the first
#           iteration this probability is exactly equal to 1); if a single number is
#           specified adaptation probability is 1 until you reach iteration equal to 
#           adapt.prob value and 0 otherwise; if a numeric vector is specified then 
#           adaptation probabilities are equal to 1 every adapt.prob[2]-th iteration 
#           unless iteration number exceeds adapt.prob[1], and 0 otherwise
#   (xi)    sample.discard: number of initial states to be discarded (length of the
#           burn-in period)
#   (xii)   thin: integer setting the thinning interval used in simulation. Only the
#           stored (thinned) points may be used for adaptation!
#   (xiii)  kern: list containing either a 'name' of a known kernel and 'par' with 
#           kernel parameters for known kernel functions (if applicable, see below); 
#           or 'fun' and 'dfun' containing functions for kernel and its derivative, 
#           respectively. Currently supported inbuilt kernels:
#           - "Gaussian": includes one scaling parameter (sigma)
#   (xiv)   log.p: logical; if TRUE algorithm works with log-densities (requires target
#           to be a log-pdf)
#   (xv)    method: string specifying the matrix decomposition used to determine the 
#           matrix root of sigma. Supported methods:
#           - and Cholesky decomposition ("chol", default); typically fastest
#           - eigenvalue decomposition ("eigen"); typically more stable
#           - singular value decomposition ("svd")
#   (xvi)   verbose: numeric value adjusting verbosity of the function. Supported
#           values:
#           - 0: no additional output is printed
#           - 1: a simple progress bar is printed
#           - 2: indicator (+) if iteration point is stored, indicator if subsample
#                has been updated (~), iteration number, current average acceptance 
#                ratio and nu parameter value, time elapsed are printed
#           - 3: same as in 2, but additionally stores calculated proposal's
#                covariance matrices
## Output
#   Following list object is returned:
#   - $x: length.out by dimension matrix containing generated MCMC output chain
#   - $accepted: vector containing the average acceptance ratio over all 
#                iterations (first element), and iterations that were stored
#                (second element)
#   - $burnin: sample.discard by dimension matrix containing the discarded
#               random walk MH chain
#   - $covMat: array containing stored proposal's covariance matrices [verbose
#              needs to be set to 3]
mcmc_kamh <- function(target, dimension=2, X0=rep(0,dimension), 
                      length.out=10000, max.subchain.length=1000, gam=0.5,
                      nu=1, nu.adapt=function(t){ifelse(t<=1000,0,1/(t-1000))}, 
                      alpha=0.234, adapt.prob=10000, sample.discard=100, thin=1,
                      kern=list(name="Gaussian",par=1), log.p=TRUE, 
                      method=c("chol","eigen", "svd"), verbose=1){
  
  ptm_startTime <- proc.time()
  
  ### Initialisations and input checks
  # Load mvtnorm for working with multivariate normal distributions
  require(mvtnorm)
  # Match target distribution
  fun_pi <- match.fun(target)
  # Set up properly adaptation probabilities function
  if(!is.numeric(adapt.prob)){
    fun_p <- match.fun(adapt.prob)
  } else {
    if(length(adapt.prob)==1){
      fun_p <- function(t){
        if(t<=adapt.prob){return(1)}
        return(0)
      }
    } else {
      if(adapt.prob[2]==1){
        adapt.prob <- adapt.prob[1]
        fun_p <- function(t){
          
          if(t<=adapt.prob){return(1)}
          return(0)}
      } else {
        fun_p <- function(t){
          if(t<=adapt.prob[1]){
            if(t%%adapt.prob[2]==1){return(1)}else{return(0)}
          } else {return(0)}  
        }       
      }
    }
  }
  if(fun_p(1)!=1){stop("Adaptation probability at the first iteration needs to be 1.")}
  # Match function for adapting nu parameter
  fun_nuAdaptScale <- match.fun(nu.adapt)
  # Set up various variables
  num_alphaStar <- alpha
  if(num_alphaStar<0|num_alphaStar>1){stop('Value of alpha needs to be between 0 and 1.')}
  num_d <- dimension
  num_nMax <- max.subchain.length
  if(num_nMax<1){Stop('Value of max.subchain.length needs to be a positive integer.')}
  num_burnin <- sample.discard + 1
  num_lengthOut <- length.out
  num_thin <- thin
  num_T <- num_lengthOut*num_thin + num_burnin
  num_gamma <- gam
  if(num_gamma==0){stop('Parameter gam cannot be zero.')}
  num_nu <- nu
  num_nuSq <- num_nu*num_nu
  mat_x <- matrix(NA, nrow=(num_burnin+num_lengthOut), ncol=num_d)
  mat_x[1,] <- X0
  num_gammaSq <- num_gamma*num_gamma
  mat_gammaSqDiag <- num_gammaSq*diag(num_d)
  num_accepted <- integer(1)
  num_acceptedAndStored <- integer(1)
  mat_z <- double(num_nMax) # Memory pre-allocation
  lst_kern <- kern
  lgc_acceptedFlag <- FALSE
  lgc_xStoredFlag <- FALSE
  lgc_xDiffersFlag <- TRUE
  lgc_nMaxNotReachedFlag <- TRUE
  # Verbose
  lgc_progressBarFlag <- (verbose==1)
  lgc_classicVerboseFlag <- (verbose>=2)
  lgc_saveCovMatFlag <- (verbose>=3)
  if(lgc_progressBarFlag){fun_progressBar <- txtProgressBar(0, num_T, initial=0, style=3)}
  if(lgc_saveCovMatFlag){
    arr_proposalCovMat <- array(NA, c(num_d, num_d, num_lengthOut))
  } else {
    arr_proposalCovMat <- "Please set verbose = 3."
  }
  
  # Determine method for generating multivariate normals
  if(num_d==1){
    # 1D case
    fun_dnorm <- function(x, mean, var){dnorm(x, mean, var*var, log=log.p)}
    fun_rnorm <- function(n, mean, var){rnorm(n, mean, var*var)}
  } else {
    # 2D or higher case
    fun_dnorm <- function(x, mean, sigmamat){dmvnorm(x, mean, sigmamat, log=log.p)}
    fun_rnorm <- function(n, mean, sigmamat){rmvnorm(n, mean, sigmamat, method=method[1])}
  }
  
  ## Detect known kernels 
  if(lst_kern$name=="Gaussian"){
    num_sigma <- lst_kern$par[1]
    num_sigmaSq <- num_sigma*num_sigma
    fun_k <- function(x,y){exp(-0.5*sum((x-y)^2)/(num_sigmaSq))}
    fun_dk <- function(x,y){exp(-0.5*sum((x-y)^2)/(num_sigmaSq))*(y-x)/num_sigmaSq}  
  } else {
    fun_k <- match.fun(lst_kern$fun)
    fun_dk <- match.fun(lst_kern$dfun)
  }

  ### Burn-in period (plus one run, giving first z)
  for(num_t in 1:(num_burnin)){
    
    ## Proposal step
    vec_xProposal <- fun_rnorm(1, mat_x[num_t,], mat_gammaSqDiag)
    
    ## Accept/Reject
    num_alpha <- ifelse(log.p==TRUE, fun_pi(vec_xProposal)-fun_pi(mat_x[num_t,]),
                        fun_pi(vec_xProposal)/fun_pi(mat_x[num_t,]))
    if(ifelse(log.p==TRUE,log(runif(1,0,1)),runif(1,0,1))<num_alpha){ 
      # Accept
      mat_x[num_t+1,] <- vec_xProposal 
    } else {
      # Reject
      mat_x[num_t+1,] <- mat_x[num_t,]
    }
    
    # Verbose
    if(lgc_progressBarFlag){setTxtProgressBar(fun_progressBar,num_t)}
    
  }
  
  ## Use the last point from burn-in phase as a first point in chain
  vec_xCurrent <- mat_x[num_burnin+1,]
  num_piAtCurrent <- fun_pi(vec_xCurrent)
  # Verbose
  if(lgc_classicVerboseFlag){
    cat("+ ",sep="")
    cat(num_t-num_burnin,". ",sep="")
    cat("(Initial step from random walk M-H burn-in phase.)")
    cat("\n",sep="")
  }

  ### Run adaptive MCMC
  if(num_lengthOut>1){for(num_t in (num_burnin+1):max((num_T-2*num_thin+1),num_burnin+1)){

    ## Subsample update
    if(runif(1,0,1)<fun_p(num_t-num_burnin)){
      # Update subsample z
      num_n <- min(ceiling((num_t-num_burnin)/thin), num_nMax)
      vec_zSampleIndices <- num_burnin + sample.int(ceiling((num_t-num_burnin)/thin), 
                                                    num_n, replace=FALSE)
      mat_z <- mat_x[vec_zSampleIndices,]
      lgc_zUpdateFlag <- TRUE
    }

    ## Calculate proposal covariance matrix
    if(num_n!=1){
      if(num_d==1){
        if(lgc_xDiffersFlag | lgc_zUpdateFlag){
          mat_M <- matrix(sapply(mat_z, function(z){2*fun_dk(vec_xCurrent, z)}), nrow=1)
        } # Optimisation: if x and z did not change do not recalculate M
      } else {
        if(lgc_xDiffersFlag | lgc_zUpdateFlag){
          mat_M <- apply(mat_z, 1, function(z){2*fun_dk(vec_xCurrent, z)})
        } # Optimisation: if x and z did not change do not recalculate M
      }
    } else {
      mat_M <- matrix(2*fun_dk(vec_xCurrent, mat_z))
    }
    
    if(lgc_nMaxNotReachedFlag){
      vec_ones <- rep(1,num_n)
      if(num_n==num_nMax){lgc_nMaxNotReachedFlag <- FALSE}
    } # Optimisation: if n does not change then do not recompute vector of ones
    mat_covAdapt <- num_nuSq*mat_M%*%t(mat_M-((mat_M%*%(vec_ones/num_n))%*%t(vec_ones)))
    
    ## Proposal step
    mat_covProposal <- mat_gammaSqDiag + mat_covAdapt
    vec_xProposal <- fun_rnorm(1, vec_xCurrent, mat_covProposal)
    if(num_t==(num_burnin+1)){
      mat_covPrevious <- mat_gammaSqDiag
    }
    
    ## Accept/Reject
    num_qzxCurrent <- fun_dnorm(vec_xCurrent, vec_xProposal, mat_covPrevious)
    num_qzxProposal <- fun_dnorm(vec_xProposal, vec_xCurrent, mat_covProposal)
    num_piAtProposal <- fun_pi(vec_xProposal)
    num_alpha <- ifelse(log.p==TRUE,
                        num_piAtProposal-num_piAtCurrent+
                          num_qzxCurrent-num_qzxProposal,
                        (num_piAtProposal/num_piAtCurrent)*
                          (num_qzxCurrent/num_qzxProposal))
    if( ifelse(log.p==TRUE, log(runif(1,0,1)), runif(1,0,1)) < num_alpha ){
      # Accept
      mat_covPrevious <- mat_covProposal
      vec_xCurrent <- as.numeric(vec_xProposal)
      num_piAtCurrent <- num_piAtProposal
      num_accepted <- num_accepted + 1
      lgc_acceptedFlag <- TRUE
      lgc_xDiffersFlag <- TRUE
    } else {
      # Do nothing for rejection - nothing changes from the previous iteration
      lgc_xDiffersFlag <- FALSE
    }   
    
    # Store the point
    if( ((num_t-num_burnin)%%num_thin==1) | (num_thin==1) ){
      mat_x[ceiling((num_t-num_burnin)/thin)+num_burnin+1,] <- vec_xCurrent
      if(lgc_saveCovMatFlag){
        arr_proposalCovMat[,,(ceiling((num_t-num_burnin)/thin))] <- mat_covProposal
      }
      lgc_xStoredFlag <- TRUE
      if(lgc_acceptedFlag){
        num_acceptedAndStored <- num_acceptedAndStored + 1
      }
    }
    
    ## Adapt nu
    num_nuAdaptScale <- fun_nuAdaptScale(num_t+1-num_burnin)
    num_nu <- exp(log(num_nu)+num_nuAdaptScale*(ifelse(lgc_acceptedFlag,1,0)-num_alphaStar))
    num_nuSq <- num_nu*num_nu
    
    ## Verbose
    if(lgc_progressBarFlag){setTxtProgressBar(fun_progressBar,num_t)}
    if(lgc_classicVerboseFlag){
      if(lgc_xStoredFlag){cat("+",sep="")}else{cat(" ",sep="")}
      if(lgc_zUpdateFlag){cat("~ ",sep="")}else{cat("  ",sep="")}
      cat(num_t-num_burnin,". ",sep="")
      cat("Acceptance ratio: ",format(round(num_accepted/(num_t-num_burnin)*100,2),
                                      nsmall=2),"% @ nu = ",
          format(round(num_nu,3),nsmall=3),". ",sep="")
      ptm_runTime <- proc.time()[1]-ptm_startTime[1] 
      ptm_estTime <- round(ptm_runTime/(num_t-num_burnin)*(num_T-num_burnin))
      cat("[Time: ",format(round(ptm_runTime,0),nsmall=0)," s]",sep="")
      cat("\n",sep="")
    }
 
    # Clean up flags
    lgc_zUpdateFlag <- FALSE
    lgc_acceptedFlag <- FALSE
    lgc_xStoredFlag <- FALSE
    
  }}#end adaptive mcmc
  
  # Verbose
  if(lgc_progressBarFlag){
    close(fun_progressBar)
    ptm_runTime <- proc.time()[1]-ptm_startTime[1] 
    cat("Done in ",ptm_runTime," s.\n",sep="")
  }
  
  ### Return results in a list
  return(list(x=mat_x[(num_burnin+1):(num_burnin+num_lengthOut),], 
              accepted=c(num_accepted, num_acceptedAndStored), 
              burnin=mat_x[1:num_burnin], covMat=arr_proposalCovMat))  
  
}