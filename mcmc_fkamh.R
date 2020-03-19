## Description
#   Fast Kernel Adaptive Metropolis-Hastings (F-KAMH) algorithm that uses random Fourier 
#   features framework to significantly improve the cost of computations by dropping
#   the dependency of calculations on the subsample size (c.f. KAMH algorithm).
## Arguments
#   (i)     target: function of the target pdf
#   (ii)    dimension: integer specifying the dimensionality of the target distribution
#   (iii)   X0: numerical vector specifying the initial state for the MCMC algorithm
#   (iv)    length.out: integer specifying output chain length, including the initial
#           state taken from the last point in the random walk M-H algorithm. Note that, 
#           this is not necessarily the number of iterations, as this specifies chain 
#           length returned after accounting for thinning.
#   (v)     gam: scaling parameter gamma, part of the proposal's covariance matrix
#   (vi)    eta: scaling parameter eta, part of the proposal's covariance matrix
#   (vii)   eta.adapt: function which takes iteration number as an argument and
#           returns the value of a scaling parameter for adapting the eta parameter
#   (viii)  alpha: target acceptance rate used in adapting the eta parameter
#   (ix)    sample.discard: number of initial states to be discarded (length of the
#           burn-in period)
#   (x)     rff.samples: number of random Fourier features (dimension D). Has to be
#           even when working with embedding = 1.
#   (xi)    thin: integer setting the thinning interval used in simulation. Only the
#           stored (thinned) points may be used for adaptation!
#   (xii)   kern: list containing either a 'name' of a known kernel, or 'fname' 
#           containing name of a known Fourier transform of the kernel, and 'par' with 
#           kernel parameters for known kernel functions (if applicable, see below);
#           or 'fun' and 'ffun' containing functions for kernel and function
#           allowing to generate i.i.d. samples from its Fourier transform, 
#           respectively. Currently supported inbuilt kernels:
#           - [name] "Gaussian": includes one scaling parameter (sigma)
#           - [fname] "Gaussian: includes one scaling parameter (sigma)
#   (xiii)  log.p: logical; if TRUE algorithm works with log-densities (requires target
#           to be a log-pdf)
#   (xiv)   method: string specifying the matrix decomposition used to determine the 
#           matrix root of sigma. Supported methods:
#           - and Cholesky decomposition ("chol", default); typically fastest
#           - eigenvalue decomposition ("eigen"); typically more stable
#           - singular value decomposition ("svd")
#   (xv)    embedding: numerical value specifying which embedding is used. Supported
#           values:
#           - 1: sine, cosine representation
#           - 2: cosine with added uniform noise b representation
#   (xvi)   verbose: numeric value adjusting verbosity of the function. Supported
#           values:
#           - 0: no additional output is printed
#           - 1: a simple progress bar is printed
#           - 2: indicator (+) if iteration point is stored, iteration number
#                current average acceptance ratio and eta parameter value, time 
#                elapsed and estimated total time required are printed
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
mcmc_fkamh <- function(target, dimension=2, X0=rep(0,dimension), 
                       length.out=10000, gam=0.5, eta=1,
                       eta.adapt=function(t){ifelse(t<=1000,0,1/(t-1000))}, 
                       alpha=0.234, sample.discard=100, rff.samples=400, 
                       thin=1, kern=list(name="Gaussian",par=1), 
                       log.p=TRUE, method=c("chol","eigen", "svd"), 
                       embedding=1, verbose=1){
  
  ptm_startTime <- proc.time()
  
  ### Initialisations and input checks
  # Load mvtnorm for working with multivariate normal distributions
  require(mvtnorm)
  # Match target distribution
  fun_pi <- match.fun(target)
  # Match function for adapting nu parameter
  fun_etaAdaptScale <- match.fun(eta.adapt)
  # Set up various variables
  num_alphaStar <- alpha
  num_embedding <- embedding
  num_d <- dimension
  num_D <- rff.samples
  if(num_embedding & num_D%%2){
    num_D <- (num_D+1)
    warning(paste("Using embedding '1' requires even number of rff.samples; ",
                  "setting dimension D to ",num_D,".",sep=""))
  }
  num_burnin <- sample.discard + 1
  num_lengthOut <- length.out
  num_thin <- thin
  num_T <- num_lengthOut*num_thin + num_burnin
  num_gamma <- gam
  if(num_gamma==0){stop('Parameter gam cannot be zero.')}
  num_eta <- eta
  num_etaSq <- num_eta*num_eta
  mat_x <- matrix(NA, nrow=(num_burnin+num_lengthOut), ncol=num_d)
  mat_x[1,] <- X0
  num_gammaSq <- num_gamma*num_gamma
  mat_gammaSqDiag <- num_gammaSq*diag(num_d)
  num_accepted <- integer(1)
  num_acceptedAndStored <- integer(1)
  lst_kern <- kern
  lgc_acceptedFlag <- FALSE
  lgc_xStoredFlag <- FALSE
  lgc_kernNotDetectedFlag <- TRUE
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

  ## Set functions in relation to the kernel
  if(!is.null(lst_kern$name)){
    # Kernel name specified; try to load known kernels
    if(lst_kern$name=="Gaussian"){
      num_sigma <- lst_kern$par[1]
      num_sigmaSq <- num_sigma*num_sigma
      fun_romega <- function(){fun_rnorm(1, rep(0,num_d), diag(num_d)/num_sigmaSq)}
      lgc_kernNotDetectedFlag <- FALSE
    } else {
      warning("Input kern$name not supported.")
    }
  } else if(lgc_kernNotDetectedFlag & (!is.null(flg_kern$fname))){
    # Fourier transform of kernel name specified; try to load known kernels
    if(lst_kern$fname=="Gaussian"){
      fun_romega <- function(){fun_rnorm(1, rep(0,num_d), diag(num_d)/lst_kern$par[1])}
      lgc_kernNotDetectedFlag <- FALSE
    } else {
      warning("Input kern$fname not supported.")
    }
  } else if(lgc_kernNotDetectedFlag) {
    fun_k <- match.fun(lst_kern$fun)
    fun_romega <- match.fun(lst_kern$ffun)
  } else {
    stop("Failed to determine the kernel function.")
  }
  
  ### Sample omega's
  mat_omegaTranspose <- apply(matrix(NA,ifelse(num_embedding==1,num_D/2,num_D),1), 1, 
                              function(x){fun_romega()})
  if(num_d==1){
    mat_omegaTranspose <- matrix(mat_omegaTranspose, nrow=1)
  }
  mat_omega <- t(mat_omegaTranspose)
  
  ### Choose appropriate embedding
  if(num_embedding==1){
    fun_phi <- function(vec_x){
      sqrt(2/num_D) * matrix(t(cbind(sin(crossprod(mat_omegaTranspose,vec_x)), 
                                     cos(crossprod(mat_omegaTranspose,vec_x)))), ncol = 1)
    }
  } else if(num_embedding==2){
    vec_b <- runif(num_D,0,2*pi)
    fun_phi <- function(vec_x){sqrt(2/num_D)*cos(crossprod(mat_omegaTranspose,vec_x)+vec_b)}
  } else {
    stop("Choice of embedding not supported.")
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
    
    ## Get embedding of current point x, vec_phiXCurrent
    vec_phiXCurrent <- fun_phi(vec_xCurrent)
    
    ## Calculate matrix of partial derivatives of phi, mat_pphi
    if(num_embedding==1){
      mat_pphi <- sqrt(2/num_D)*matrix(t(
        cbind(sweep(mat_omega,MARGIN=1,cos(rowSums(sweep(mat_omega,MARGIN=2,
                                                         vec_xCurrent,"*"))), "*"),
              -sweep(mat_omega,MARGIN=1,sin(rowSums(sweep(mat_omega,MARGIN=2,
                                                          vec_xCurrent,"*"))), "*")
              )
        ), ncol=num_D, byrow=FALSE)
    } else {
      mat_pphi <- -sqrt(2/num_D)*sweep(mat_omega, 
                                       MARGIN=1, 
                                       sin(rowSums(sweep(mat_omega,MARGIN=2,
                                                         vec_xCurrent,"*")) + vec_b), "*")
    }
    
    ## Calculate matrix C using rank-one update
    if(num_t==(num_burnin+1)){
      vec_muCurrent <- vec_phiXCurrent
      mat_M <- matrix(0, num_D, num_D)
      mat_C <- matrix(0, num_D, num_D)
      mat_covPrevious <- mat_gammaSqDiag
    } else {
      vec_muPrevious <- vec_muCurrent
      vec_muCurrent <- ((num_t-num_burnin)/((num_t-num_burnin)+1))*vec_muPrevious + 
        vec_phiXCurrent/((num_t-num_burnin)+1)
      mat_M <- mat_M + (vec_phiXCurrent - vec_muPrevious) %*% t(vec_phiXCurrent - 
                                                                  vec_muCurrent)
      mat_C <- mat_M/((num_t-num_burnin)+1)
    }
    
    ## Proposal step
    if(num_embedding==1){
      mat_covAdapt <- num_etaSq*mat_pphi%*%mat_C%*%t(mat_pphi)
    } else {
      mat_covAdapt <- num_etaSq*t(mat_pphi)%*%mat_C%*%mat_pphi
    }
    mat_covProposal <- mat_gammaSqDiag + mat_covAdapt
    vec_xProposal <- fun_rnorm(1, vec_xCurrent, mat_covProposal)

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
    }
    # Do nothing for rejection - nothing changes from the previous iteration
    
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
    
    ## Adapt eta
    num_etaAdaptScale <- fun_etaAdaptScale(num_t+1-num_burnin)
    num_eta <- exp(log(num_eta)+num_etaAdaptScale*(ifelse(lgc_acceptedFlag,1,0)-
                                                     num_alphaStar))
    num_etaSq <- num_eta*num_eta
    
    ## Verbose
    if(lgc_progressBarFlag){setTxtProgressBar(fun_progressBar,num_t)}
    if(lgc_classicVerboseFlag){
      if(lgc_xStoredFlag){cat("+ ",sep="")}else{cat("  ",sep="")}
      cat(num_t-num_burnin,". ",sep="")
      cat("Acceptance ratio: ",format(round(num_accepted/(num_t-num_burnin)*100,2),
                                      nsmall=2),"% @ eta = ",
          format(round(num_eta,3),nsmall=3),". ",sep="")
      ptm_runTime <- proc.time()[1]-ptm_startTime[1] 
      ptm_estTime <- round(ptm_runTime/(num_t-num_burnin)*(num_T-num_burnin))
      cat("[Time: ",format(round(ptm_runTime,0),nsmall=0)," / ",
          format(round(ptm_estTime,0),nsmall=0)," s]",sep="")
      cat("\n",sep="")
    }

    # Clean up flags
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