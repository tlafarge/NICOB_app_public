######################################################################
##
## FILE        : DoEUnilateralBayes.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Aug-02

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
##         mcmc = MCMC output from bayesGelman is used, if
##                supplied, only when LOO=FALSE (optional)
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

DoEUnilateralBayes = function (x, u, nu, lab, LOO, mcmc,
                               ni, nb, nt,coverageProb,##AK change## removed default values and added variable coverage
                               UItauPriorScale, UIsigmaPriorScale)  ##AK change## no longer calculate new hyper-parameters for LOO
{
  source("bayesGelman.R") ##AK change##

  n = length(x)

  ## check convergence of LOO runs
  bayesRunsConv=c()

  consensusValueName="mu" ##AK change##
  StudyEffectName="xi" ##AK change##

  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}

  # Check if any lab was exempt from prior Mcmc computation
  sanitize = !startsWith( lab,"-")

  #Number of labs taken into account for the consensus
  nI = sum(sanitize)

  if (LOO) {
    nc = length(mcmc)
    ## DoEs based on leave-one-out estimates
    cat("## DoEUnilateralBayes (Leave-One-Out Version) ---\n")

    ##Add Warning if LOO is one with only one lab used for the consensusV
    if (nI != n && nI<= 1){
      testWarnLab="WARNING: Not enough labs were including in the consensus for DoE with LOO to be meaningful<br/>"
    }

    # ##### Code below is from Antonio's 2016-Nov-02 modification of the DoE code
    # ##### I don't use sigmaMCMC in the LOO version of the DoE, per "Questions and Answers" emails (August 2016).
    # ## We begin by populating sigmaMCMC that will be used
    # ## below (to define gammaj) using the results of a
    # ## Bayes-Gelman model fitted to all of the data
    #
    # # if (is.null(mcmc)) {mcmc = bayesGelman(x=x, u=u, nu=nu, ni=ni, nb=nb, nt=nt)$mcmc} ##AK change## this option unnecessary for NICOB
    #
    #
    # if (is.null(nu)) {
    #   ## When nu is NULL we take sigma=u, as if all elements
    #   ## of u were based on infinitely many degrees of
    #   ## freedom, and the sigma will not have been estimated
    #   K = sum(sapply(mcmc, nrow))
    #   sigmaMCMC =  t(array(rep(u, K), dim=c(n,K)))
    # } else {
    #   ## When nu is not NULL, the sigma will have
    #   ## been estimated and there are samples from
    #   ## their posterior distributions
    #   sigmaMCMC = NULL
    #   # sigmaNAMEs = paste("sigma.", lab, sep="") ##AP labels##
    #   sigmaNAMEs = paste("sigma[", 1:n,"]", sep="") ##AK change##
    #
    #   for (jc in 1:nc)
    #   {
    #     is = match(sigmaNAMEs, dimnames(mcmc[[jc]])[[2]])
    #     sigmaMCMC = rbind(sigmaMCMC, mcmc[[jc]][,is])
    #   }
    # }

    ## For each j=1,...,n, BG[[j]] is a matrix with two
    ## columns and K rows: the columns have samples of the
    ## posterior distributions of mu and tau when the
    ## measurement results from lab j were st aside
    BG = list()
    for (j in 1:n)
    {
      cat(j, "of", n, "\n")
      zj = bayesGelman(x=x[-j], u=u[-j], nu=nu[-j], ##AK change## don't use lab labels in NICOB bayesGelman
                       tauPriorScale=UItauPriorScale,#mad(x[-j]),
                       sigmaPriorScale=UIsigmaPriorScale,#median(u[-j]),
                       nc=1,
                       ni=ni, nb=nb, nt=nt)
      bayesRunsConv=c(bayesRunsConv,(zj$warn!=""))

      BG[[j]] = rbind(zj$mcmc[[1]][, c(consensusValueName, "tau")]) ##AK change## only one chain now
                      # zj$mcmc[[2]][, c(consensusValueName, "tau")])
    }

    K = nrow(BG[[1]])
    D = array(dim=c(K, n))
    DoE.x = DoE.U = DoE.Lwr = DoE.Upr = numeric(n)
    for (j in 1:n)
    {
      muMCMC = BG[[j]][,consensusValueName] ##AK change##
      tauMCMC = BG[[j]][,"tau"]

      ######## Used below to generate ej, different from 2016-Nov-02 update
      sigmaj = sqrt(tauMCMC^2 + u[j]^2)
      if (is.null(nu[j]) || (nu[j]==Inf)) {
      # if (is.null(nu[j]) | (isTRUE(nu[j]==Inf))) {
        ej = rnorm(K, mean=0, sd=sigmaj)
      } else {
        ## AP-AK: Addressing possibility of no
        ## more than 2 degrees of freedom
        phi = if (nu[j] <= 2) {## other options: ## sqrt(3) ##qt(.975,df=nu[j])/qt(.975,df=Inf)
          a = 0.4668994; b= -0.3998882
          1.5*((1 - (3/4)*(a - 4*b *
                             (nu[j]^(-3/4)-1)/3))^(-4/3))
        } else {sqrt(nu[j]/(nu[j]-2))}
        # ej = sigmaj*u[j]*rt(K, df=nu[j])/phi ##AK change## typo
        ej = sigmaj*rt(K, df=nu[j])/phi

      }
      ########

      # ######## Generation of ej from 2016-Nov-02 update
      # gammaj = sqrt(tauMCMC^2 + sigmaMCMC[,j]^2)
      # ej = rnorm(K, mean=0, sd=gammaj)
      # ########

      DoE.x[j] = x[j] - mean(muMCMC)
      D[,j] = x[j] - muMCMC + ej
      DoE.U[j] =
        symmetricalBootstrapCI(D[,j], mean(D[,j]), coverageProb) ##AK change##
      DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)#0.025) ##AK change##
      DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)#0.975) ##AK change##
    }

  } else {

    ## DoEs computed according to MRA


    nc = length(mcmc)
    if (is.null(nu)) {
      ## When nu is NULL we take sigma=u, as if all elements
      ## of u were based on infinitely many degrees of
      ## freedom, and the sigma will not have been estimated
      muMCMC = tauMCMC = NULL
      for (jc in 1:nc)
      {
        muMCMC = c(muMCMC, mcmc[[jc]][,consensusValueName]) ##AK change##
        tauMCMC = c(tauMCMC, mcmc[[jc]][,"tau"])
      }

      K = length(muMCMC)

      sigmaMCMC =  t(array(rep(u, K), dim=c(n,K)))
      } else {
        ## When nu is not NULL, the sigma will have
        ## been estimated and there are samples from
        ## their posterior distributions
        muMCMC = thetaMCMC = tauMCMC = sigmaMCMC = NULL

        # thetaNAMEs = paste("theta.", lab, sep="")  ##AP labels##
        # sigmaNAMEs = paste("sigma.", lab, sep="")  ##AP labels##
        thetaNAMEs = paste(StudyEffectName, "[", 1:n,"]", sep="") ##AK change##
        sigmaNAMEs = paste("sigma[", 1:n,"]", sep="") ##AK change##


        for (jc in 1:nc)
        {
          muMCMC = c(muMCMC, mcmc[[jc]][,consensusValueName]) ##AK change##
          it = match(thetaNAMEs, dimnames(mcmc[[jc]])[[2]])
          thetaMCMC = rbind(thetaMCMC, mcmc[[jc]][,it])
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
    ## obtain (M,T,S), and then making a drawing Y from a
    ## Gaussian distribution with mean M and standard
    ## deviation sqrt(T^2+S^2). The corresponding value of D
    ## computed below is Y-M.  The unilateral DoE for lab j is
    ## x[j]-muHAT, where muHAT is the mean of the posterior
    ## distribution of mu. The U for this DoE is U(D)
    D = array(dim=c(K, n))
    DoE.x = DoE.U = DoE.Lwr = DoE.Upr = numeric(n)


      #We compute the DoE for the lab included in the consensus first
      for (i in 1:nI)
      {

        j = (1:n)[sanitize][i]
        D[,j] = x[j] - muMCMC +
          rnorm(K, mean=0, sd=sqrt(tauMCMC^2 + sigmaMCMC[,i]^2))
        DoE.x[j] = x[j] - muHAT
        DoE.U[j] =
          symmetricalBootstrapCI(D[,j], mean(D[,j]), coverageProb) ##AK change##
        DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)#0.025) ##AK change##
        DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)#0.975) ##AK change##
      }

      #We compute the DoE for the lab exclude from the consensus
      #We use a formula in case nu is specified that requires 3 or more labs to be included in the consensus

      if(nI<n)
      {
        if (is.null(nu))
          nu=rep(Inf,n)

        for (i in 1:(n-nI))
        {

          j = (1:n)[!sanitize][i]
          nutau= nI-1
          nudj=( u[j]^2/(nu[j]+1) + tauMCMC^2/(nutau+1) )^2 / (  (u[j]^2/(nu[j]+1))^2/nu[j]  + (tauMCMC^2/(nutau+1))^2/nutau   )
          D[,j] = x[j] - muMCMC +      sqrt(tauMCMC^2 +u[j]^2)/sqrt(nudj/(nudj-2))*rt(K, nudj)
          DoE.x[j] = x[j] - muHAT
          DoE.U[j] =
            symmetricalBootstrapCI(D[,j], mean(D[,j]), coverageProb) ##AK change##
          DoE.Lwr[j] = quantile(D[,j], probs=(1-coverageProb)/2)#0.025) ##AK change##
          DoE.Upr[j] = quantile(D[,j], probs=(1+coverageProb)/2)#0.975) ##AK change##


        }

      }
    }


  results = data.frame(Lab=lab, DoE.x=DoE.x, DoE.U95=DoE.U,
                       DoE.Lwr=DoE.Lwr, DoE.Upr=DoE.Upr)


  ### Send warning if MCMC didn't converge
  if (sum(bayesRunsConv)>0) {
    testWarn=paste("WARNING: MCMC may not have reached equilibrium<br/>",
                   "WARNING: Results are unreliable<br/>",
                   paste("SUGGESTION: Re-run bayesGelman with ",
                         "ni = ", 10*ni, ", nb = ", 2*nb,
                         ", and nt = ", ceiling(0.001*ni), "<br/>", sep=""),sep=" ")

  }else{
    testWarn=""
  }
  if (exists("testWarnLab"))
  {
    testWarn = paste(testWarnLab, testWarn)
  }

  return(invisible(list(D=D, DoE=results,DoEwarn=testWarn)))
}

## ######################################################################
## ##
## ## EXAMPLE
##
## source("~/NIST/TALES/Analysis/symmetricalBootstrapCI.R")
## source("~/NIST/TALES/Analysis/bayesGelman.R")
##
## PCB28 = data.frame(lab=c("IRMM", "KRISS", "NARL", "NIST", "NMIJ", "NRC"),
##     x=c(34.3, 32.9, 34.53, 32.42, 31.9, 35.8),
##     u=c(1.03, 0.69, 0.83, 0.29, 0.4, 0.38),
##     nu=c(60, 4, 18, 2, 13, 60))
##
## PCB28.Bayes = bayesGelman(x=PCB28$x, u=PCB28$u, nu=PCB28$nu, lab=PCB28$lab)
## PCB28.Bayes.DoEUni.MRA =
##   DoEUnilateralBayes(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##                      lab=PCB28$lab, mcmc=PCB28.Bayes$mcmc, LOO=FALSE)
##
## ##   lab      DoE.x  DoE.U95 DoE.Lwr95 DoE.Upr95
## ##  IRMM  0.6942624 4.539187 -3.837794  5.245206
## ## KRISS -0.7057376 4.313178 -5.004783  3.620437
## ##  NARL  0.9242624 4.287123 -3.448236  5.128103
## ##  NIST -1.1857376 4.132620 -5.329335  2.939262
## ##  NMIJ -1.7057376 4.098493 -5.730298  2.480700
## ##   NRC  2.1942624 4.074045 -1.935093  6.214643
##
## PCB28.DoEUni.LOO = DoEUnilateralBayes(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##      lab=PCB28$lab, LOO=TRUE)
##
## ##   lab      DoE.x  DoE.U95  DoE.Lwr95 DoE.Upr95
## ##  IRMM  0.8017541 4.929384 -4.1309793  5.723298
## ## KRISS -0.8528847 5.092199 -5.9882476  4.155381
## ##  NARL  1.0907891 4.801484 -3.7783236  5.849077
## ##  NIST -1.4618595 4.662897 -6.3052350  3.068161
## ##  NMIJ -2.1237910 4.144308 -6.1387435  2.135147
## ##   NRC  2.8500963 2.923678 -0.1168719  5.711596

######################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################
