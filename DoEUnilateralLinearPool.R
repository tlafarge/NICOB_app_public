######################################################################
##
## FILE        : DoEUnilateralLinearPool.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Jun-05

## INPUT : x.All   = Numeric vector with values measured by the labs
##         u.All   = Numeric vector with standard uncertainties
##               associated with measured values
##         nu.All  = Numeric vector with (positive) numbers of
##               degrees of freedom that the elements of u are
##               based on; if NULL, then the elements of u are
##               assumed to be based on infinitely many degrees
##               of freedom; any elements of nu set to Inf are
##               replaced by 1000
##         lab.All = Character vector with laboratory labels
##         weights = weights for included labs
##         K   = Number of bootstrap replicates (optional)
##         LOO = TRUE/FALSE indicating whether to base the DoE on
##               differences to leave-one-out estimates of the
##               consensus value, or on differences to the
##               consensus value that has been derived from all
##               measurement results (optional)
##         linearOpRes = The results of the linear pooling function

## OUTPUT: List with two elements: (D) Kxn array whose jth column
## has K bootstrap replicates of the unilateral DoE for lab j;
## (DoE) data frame with values and 95 % expanded uncertainties
## for unilateral degrees of equivalence, and 2.5th and 97.5th
## percentiles of the bootstrap distributions of these values

######################################################################

DoEUnilateralLinearPool = function (x.All, u.All, nu.All, lab.All, weights,
                                    K, LOO, coverageProb, linearOpRes) 
{

  # source("linearOP.R") #already ran in server.R
  testWarn = ""
  
  n.All = length(x.All)
  if (is.null(lab.All)) {lab.All=paste("L", 1:n.All, sep="")}

  # Check if any lab was exempt from prior computation
  sanitize = !startsWith( lab.All,"-")
  
  ## Analysis for included labs 
  x=x.All[sanitize]
  u=u.All[sanitize]
  nu=nu.All[sanitize]
  lab=lab.All[sanitize]

  #Number of labs taken into account for the consensus
  nI = sum(sanitize)
  if (nI != n.All && nI<= 2 && LOO){
    testWarn="WARNING: Not enough labs were including in the consensus for DoE with LOO to be meaningful<br/>"
  }

  z.POOL = linearOpRes
  mu = mean(z.POOL)
  DoE.x = DoE.Lwr = DoE.Upr = DoE.U = numeric(nI)
  names(DoE.x) = names(DoE.U) = lab
  names(DoE.Lwr) = names(DoE.Upr) = lab
  D = array(dim=c(K,nI))
  dimnames(D)[[2]] = lab

  if (LOO) {
    ## DoEs based on leave-one-out estimates
    cat("## DoEUnilateralLinearPool (LOO) --------------\n")

    # # Check if any lab was exempt from prior Mcmc computation
    # sanitize = !startsWith( lab,"-")
    # 
    # #Number of labs taken into account for the consensus
    # nI = sum(sanitize)
    # if (nI != n && nI<= 1){
    #   testWarn="WARNING: Not enough labs were including in the consensus for DoE with LOO to be meaningful<br/>"
    # }

    for (j in 1:nI)
    {

      zj = linearOP(x=x[-j], u=u[-j], nu=nu[-j], weights=weights[-j], K)  
      DoE.x[j] = x[j] - mean(zj)
      # ## OPTION A
      # ## ej plays the role of a sample drawn from
      # ## the distribution of laboratory effects in
      # ## a random effects model
      # ej = sample(zj-mean(zj), size=K, replace=TRUE)
      ## OPTION B
      ## A more conservative alternative would be
      ## ej = rnorm(K, mean=0, sd=u[j]), or a
      ## corresponding sample drawn from a
      ## re-scaled Student's t(nu[j])
      sigmaj = u[j]
      if (is.null(nu[j]) || (nu[j]==Inf)) {
        ej = rnorm(K, mean=0, sd=sigmaj)
      } else {
        ## AP-AK: Addressing possibility of no
        ## more than 2 degrees of freedom
        phi = if (nu[j] <= 2) {## other options: ## sqrt(3) ##qt(.975,df=nu[j])/qt(.975,df=Inf)
          a = 0.4668994; b= -0.3998882
          1.5*((1 - (3/4)*(a - 4*b *
                             (nu[j]^(-3/4)-1)/3))^(-4/3))
        } else {sqrt(nu[j]/(nu[j]-2))}
        ej = sigmaj*rt(K, df=nu[j])/phi

      }

      D[,j] = x[j] + ej - zj
      DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)#0.025)
      DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)#0.975)
      DoE.U[j] = symmetricalBootstrapCI(D[,j], DoE.x[j], coverageProb)
    }

  } else {

    ## DoEs computed according to MRA
    cat("## DoEUnilateralLinearPool (MRA) --------------\n")
    for (j in 1:nI)
    {

      if (is.null(nu[j]) || (nu[j]==Inf)) {
        ej = rnorm(K, mean=0, sd=u[j])
      }else{

        ## AP-AK: Addressing possibility of no
        ## more than 2 degrees of freedom
        phi = if (nu[j] <= 2) {## other options: ## sqrt(3) ##qt(.975,df=nu[j])/qt(.975,df=Inf)
          a = 0.4668994; b= -0.3998882
          1.5*((1 - (3/4)*(a - 4*b *
                             (nu[j]^(-3/4)-1)/3))^(-4/3))
        } else {sqrt(nu[j]/(nu[j]-2))}

        factorVar = u[j]/phi

        ej = rt(K, df=nu[j]) * factorVar

      }

      DoE.x[j] = x[j] - mu
      D[,j] = x[j] + ej - z.POOL
      DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)#0.025)
      DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)#0.975)
      DoE.U[j] = symmetricalBootstrapCI(D[,j], DoE.x[j], coverageProb)
    }

  } ### end of included lab analysis
  
  

  ### Analysis for excluded labs
  
  if(n.All>nI){
    x.excluded=x.All[!sanitize]
    u.excluded=u.All[!sanitize]
    nu.excluded=nu.All[!sanitize]
    lab.excluded=lab.All[!sanitize]
    
    DoE.x.excluded = DoE.Lwr.excluded = DoE.Upr.excluded = DoE.U.excluded = numeric(n.All-nI)
    names(DoE.x.excluded) = names(DoE.U.excluded) = names(DoE.Lwr.excluded) = names(DoE.Upr.excluded) =lab.excluded
    
    DoE.U.excluded = rep(NA,n.All-nI)
    D.excluded = array(NA,dim=c(K,n.All-nI)) 
    
    dimnames(D.excluded)[[2]] = as.list(lab.excluded)
    
    for(j in 1:(n.All-nI)){
      
      #### AP CHECK: using same code from MRA, looks reasonable 
      if (is.null(nu.excluded[j]) || (nu.excluded[j]==Inf)) {
        ej = rnorm(K, mean=0, sd=u.excluded[j])
      }else{
        
        phi = if (nu.excluded[j] <= 2) {
          a = 0.4668994; b= -0.3998882
          1.5*((1 - (3/4)*(a - 4*b *
                             (nu.excluded[j]^(-3/4)-1)/3))^(-4/3))
        } else {sqrt(nu.excluded[j]/(nu.excluded[j]-2))}
        
        factorVar = u.excluded[j]/phi
        
        ej = rt(K, df=nu.excluded[j]) * factorVar
        
      }
      
      DoE.x.excluded[j] = x.excluded[j] - mu
      D.excluded[,j] = x.excluded[j] + ej - z.POOL
      DoE.Lwr.excluded[j] = quantile(D.excluded[,j], probs=(1-coverageProb)/2)
      DoE.Upr.excluded[j] = quantile(D.excluded[,j], probs=(1+coverageProb)/2)
      DoE.U.excluded[j] = symmetricalBootstrapCI(D.excluded[,j], DoE.x.excluded[j], coverageProb)
      
    }  

    DoE.x=c(DoE.x,DoE.x.excluded)
    DoE.U=c(DoE.U,DoE.U.excluded)
    DoE.Lwr=c(DoE.Lwr,DoE.Lwr.excluded) 
    DoE.Upr=c(DoE.Upr,DoE.Upr.excluded)
    
    lab.outlabel=c(lab,lab.excluded)
    
    D=cbind(D,D.excluded)
    
  }else{
    lab.outlabel=lab.All
  }
  
  ### ### ### ### 
  
  results = data.frame(Lab=lab.outlabel, DoE.x=DoE.x, DoE.U95=DoE.U,
                       DoE.Lwr=DoE.Lwr, DoE.Upr=DoE.Upr)
  row.names(results) = NULL
  # if (printResults) {print(results)}
  return(invisible(list(D=D, DoE=results,DoEwarn=testWarn)))
}

## ######################################################################
## ##
## ## EXAMPLE
##
## source("~/NIST/TALES/Analysis/linearPool.R")
## source("~/NIST/TALES/Analysis/symmetricalBootstrapCI.R")
##
## PCB28 = data.frame(lab=c("IRMM", "KRISS", "NARL", "NIST", "NMIJ", "NRC"),
##     x=c(34.3, 32.9, 34.53, 32.42, 31.9, 35.8),
##     u=c(1.03, 0.69, 0.83, 0.29, 0.4, 0.38),
##     nu=c(60, 4, 18, 2, 13, 60))
##
## DoE.LinPool.MRA = DoEUnilateralLinearPool(x=PCB28$x, u=PCB28$u,
##     nu=PCB28$nu, lab=PCB28$lab, LOO=FALSE)
##
## ##   lab      DoE.x  DoE.U95  DoE.Lwr95 DoE.Upr95
## ##  IRMM  0.6600223 3.398511 -2.8854658 3.8851724
## ## KRISS -0.7399777 2.930887 -3.8154223 1.9742119
## ##  NARL  0.8900223 3.137559 -2.4044592 3.8301882
## ##  NIST -1.2199777 2.640942 -4.0040545 1.1582162
## ##  NMIJ -1.7399777 2.682553 -4.5784838 0.6842007
## ##   NRC  2.1600223 2.670544 -0.6704649 4.5642498
##
## DoE.LinPool.MRA = DoEUnilateralLinearPool(x=PCB28$x, u=PCB28$u,
##     nu=PCB28$nu, lab=PCB28$lab, LOO=TRUE)
##
## ## OPTION A
## ##   lab      DoE.x  DoE.U95  DoE.Lwr95 DoE.Upr95
## ##  IRMM  0.7930840 4.059459 -3.2697537  4.848974
## ## KRISS -0.8891831 4.083265 -4.9690554  3.196902
## ##  NARL  1.0663682 4.078318 -3.0093121  5.147333
## ##  NIST -1.4673717 4.055848 -5.5247796  2.586661
## ##  NMIJ -2.0904674 3.658127 -5.7499717  1.566514
## ##   NRC  2.5918255 3.445461 -0.8537822  6.037205
##
## ## OPTION B
## ##   lab      DoE.x  DoE.U95  DoE.Lwr95 DoE.Upr95
## ##  IRMM  0.7908939 3.466462 -2.8946847 3.9740915
## ## KRISS -0.8905151 2.941449 -3.8977047 1.9765109
## ##  NARL  1.0663655 3.240230 -2.4209599 3.9220865
## ##  NIST -1.4651270 2.613303 -4.0639983 1.1638013
## ##  NMIJ -2.0899122 2.464111 -4.6470850 0.2407082
## ##   NRC  2.5901485 2.449019 -0.1567995 4.6254021

######################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################
