######################################################################
##
## FILE        : DoEUnilateralBayesLaplace.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2018-May-29

## INPUT : x    = Numeric vector with values measured by the labs
##         u    = Numeric vector with standard uncertainties
##                associated with measured values
##         nu   = Numeric vector with (positive) numbers of
##                degrees of freedom that the elements of u are
##                based on; if NULL, then the elements of u are
##                assumed to be based on infinitely many degrees
##                of freedom (optional)
##         lab  = Character vector with laboratory labels (optional)
##         LOO  = TRUE/FALSE indicating whether to base the DoE on
##                differences to leave-one-out estimates of the
##                consensus value, or on differences to the
##                consensus value that has been derived from all
##                measurement results (optional)
##         mcmc = MCMC output from bayesLaplace is used only when LOO=FALSE
##         ni   = Number of total MCMC draws for each chain
##         nb   = Number of MCMC draws to discard at the beginning
##                of each chain
##         nt   = Thinning rate (positive integer) for each chain

## OUTPUT: List with two elements: (D) Kxn array whose jth column
## has K bootstrap replicates of the unilateral DoE for lab j;
## (DoE) data frame with values and 95 % expanded uncertainties
## for unilateral degrees of equivalence, and 2.5th and 97.5th
## percentiles of the bootstrap distributions of these values

######################################################################

DoEUnilateralBayesLaplace = function (x, u, nu, lab, LOO, mcmc,
                               ni, nb, nt,coverageProb,
                               UItauPriorScale, UIsigmaPriorScale) 
{
  source("bayesLaplace.R") 
  
  n = length(x)
  
  ## check convergence of LOO runs
  bayesRunsConv=c() 
  
  consensusValueName="mu" 
  StudyEffectName="xi" 
  
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}
  
  if (LOO) {
    ## DoEs based on leave-one-out estimates
    cat("## DoEUnilateralBayesLaplace (Leave-One-Out Version) ---\n")
    
    ## We begin by populating sigmaMCMC that will be used
    ## below (to define gammaj) using the results of a
    ## Bayes-Laplace model fitted to all of the data
    nc = length(mcmc)
    if (is.null(nu)) {
      ## When nu is NULL we take sigma=u, as if all elements
      ## of u were based on infinitely many degrees of
      ## freedom, and the sigma will not have been estimated
      K = sum(sapply(mcmc, nrow))
      sigmaMCMC =  t(array(rep(u, K), dim=c(n,K)))
    } else {
      ## When nu is not NULL, the sigma will have
      ## been estimated and there are samples from
      ## their posterior distributions
      sigmaMCMC = NULL
      sigmaNAMEs = paste("sigma[", 1:n,"]", sep="")
      
      for (jc in 1:nc)
      {
        is = match(sigmaNAMEs, dimnames(mcmc[[jc]])[[2]])
        sigmaMCMC = rbind(sigmaMCMC, mcmc[[jc]][,is])
      }
    }

    ## For each j=1,...,n, BL[[j]] is a matrix with two
    ## columns and K rows: the columns have samples of the
    ## posterior distributions of mu and tau when the
    ## measurement results from lab j were st aside
    BL = list()
    for (j in 1:n)
    {
      cat(j, "of", n, "\n")
      zj = bayesLaplace(x=x[-j], u=u[-j], nu=nu[-j], 
                       tauPriorScale=UItauPriorScale,#mad(x[-j]),
                       sigmaPriorScale=UIsigmaPriorScale,#median(u[-j]), 
                       nc=1,
                       ni=ni, nb=nb, nt=nt)
      bayesRunsConv=c(bayesRunsConv,(zj$warn!=""))
      
      BL[[j]] = rbind(zj$mcmc[[1]][, c(consensusValueName, "tau")]) 
    }
    
    ## NOTE: This K applies to the case being entertained,
    ## where the elements of BL are based on the MCMC sample
    ## values from 2 chains exactly (nc=2) -- refer to the
    ## previous NOTE
    K = nrow(BL[[1]])
    D = array(dim=c(K, n))
    DoE.x = DoE.U = DoE.Lwr = DoE.Upr = numeric(n)
    for (j in 1:n)
    {
      muMCMC = BL[[j]][,consensusValueName] 
      tauMCMC = BL[[j]][,"tau"]

      e1j = rlaplace(K, m=0, s = tauMCMC/sqrt(2))
      e2j = rnorm(K, mean = 0, sd = sigmaMCMC)
      ej=e1j+e2j
      
      DoE.x[j] = x[j] - mean(muMCMC)
      D[,j] = x[j] - muMCMC + ej
      DoE.U[j] =
        symmetricalBootstrapCI(D[,j], mean(D[,j]), coverageProb) 
      DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)
      DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)
    }
    
  } else {
    
    ## DoEs computed according to MRA

    nc = length(mcmc)
    if (is.null(nu)) {
      ## When nu is NULL we take sigma=u, as if all elements
      ## of u were based on infinitely many degrees of
      ## freedom, and the sigma will not have been estimated
      K = sum(sapply(mcmc, nrow))
      sigmaMCMC =  t(array(rep(u, K), dim=c(n,K)))
      
      muMCMC = tauMCMC = NULL
      for (jc in 1:nc)
      {
        muMCMC = c(muMCMC, mcmc[[jc]][,consensusValueName]) 
        tauMCMC = c(tauMCMC, mcmc[[jc]][,"tau"])
      } 
      
    } else {
      ## When nu is not NULL, the sigma will have
      ## been estimated and there are samples from
      ## their posterior distributions
      
      muMCMC = tauMCMC = sigmaMCMC = NULL
      
      sigmaNAMEs = paste("sigma[", 1:n,"]", sep="")  
      
      for (jc in 1:nc)
      {
        muMCMC = c(muMCMC, mcmc[[jc]][,consensusValueName]) 
        tauMCMC = c(tauMCMC, mcmc[[jc]][,"tau"])
        is = match(sigmaNAMEs, dimnames(mcmc[[jc]])[[2]])
        sigmaMCMC = rbind(sigmaMCMC, mcmc[[jc]][,is])
      }
    }
    K = length(muMCMC)
    muHAT = mean(muMCMC)
    ## A drawing from the predictive distribution for lab j is
    ## obtained by first making a drawing from the joint
    ## posterior distribution of (mu, tau, sigma[j]), to
    ## obtain (M,T,S). The corresponding value of D
    ## computed below is Y-M+e1j+e2j, where e1j and e2j are random 
    ## draws from a Laplace and Gaussian distribution with standard 
    ## deviations corresponding to T and S.  The unilateral DoE for lab j is
    ## x[j]-muHAT, where muHAT is the mean of the posterior
    ## distribution of mu. The U for this DoE is U(D)
    D = array(dim=c(K, n))
    
    DoE.x = DoE.U = DoE.Lwr = DoE.Upr = numeric(n)
    for (j in 1:n)
    {
      D[,j] = x[j] - muMCMC + 
        rlaplace(K,m = 0, s = tauMCMC/sqrt(2)) +
        rnorm(K, mean=0, sd=sigmaMCMC[,j])
      DoE.x[j] = x[j] - muHAT
      DoE.U[j] =
        symmetricalBootstrapCI(D[,j], mean(D[,j]), coverageProb) 
      DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)
      DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)
    }
  }
  
  results = data.frame(Lab=lab, DoE.x=DoE.x, DoE.U95=DoE.U,
                       DoE.Lwr=DoE.Lwr, DoE.Upr=DoE.Upr)
  
  
  ### Send warning if MCMC didn't converge
  if (sum(bayesRunsConv)>0) {
    testWarn=paste("WARNING: MCMC may not have reached equilibrium<br/>",
                   "WARNING: Results are unreliable<br/>",
                   paste("SUGGESTION: Re-run with ",
                         "ni = ", 10*ni, ", nb = ", 2*nb,
                         ", and nt = ", ceiling(0.001*ni), "<br/>", sep=""),sep=" ") 
    
  }else{
    testWarn=""
  }
  
  return(invisible(list(D=D, DoE=results,DoEwarn=testWarn)))
}

