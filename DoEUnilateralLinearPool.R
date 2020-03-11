######################################################################
##
## FILE        : DoEUnilateralLinearPool.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Jun-05

## INPUT : x   = Numeric vector with values measured by the labs
##         u   = Numeric vector with standard uncertainties
##               associated with measured values
##         nu  = Numeric vector with (positive) numbers of
##               degrees of freedom that the elements of u are
##               based on; if NULL, then the elements of u are
##               assumed to be based on infinitely many degrees
##               of freedom; any elements of nu set to Inf are
##               replaced by 1000
##         lab = Character vector with laboratory labels
##         K   = Number of bootstrap replicates (optional)
##         LOO = TRUE/FALSE indicating whether to base the DoE on
##               differences to leave-one-out estimates of the
##               consensus value, or on differences to the
##               consensus value that has been derived from all
##               measurement results (optional)
##         linearOpRes = The results of the linear pooling fonction
##         printResults = TRUE/FALSE indicating whether to print
##               table of results (optional)

## OUTPUT: List with two elements: (D) Kxn array whose jth column
## has K bootstrap replicates of the unilateral DoE for lab j;
## (DoE) data frame with values and 95 % expanded uncertainties
## for unilateral degrees of equivalence, and 2.5th and 97.5th
## percentiles of the bootstrap distributions of these values

######################################################################

DoEUnilateralLinearPool = function (x, u, nu, lab, weights,
                                    K, LOO,coverageProb, linearOpRes) ##AK change## removed default LOO setting and added variable coverage
{

  source("linearOP.R") ##AK change##
  testWarn = ""
  n = length(x)
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}
  z.POOL = linearOpRes
  mu = mean(z.POOL)
  DoE.x = DoE.Lwr = DoE.Upr = DoE.U = numeric(n)
  names(DoE.x) = names(DoE.U) = lab
  names(DoE.Lwr) = names(DoE.Upr) = lab
  D = array(dim=c(K,n))
  dimnames(D)[[2]] = lab

  if (LOO) {
    ## DoEs based on leave-one-out estimates
    cat("## DoEUnilateralLinearPool (LOO) --------------\n")

    # Check if any lab was exempt from prior Mcmc computation
    sanitize = !startsWith( lab,"-")

    #Number of labs taken into account for the consensus
    nI = sum(sanitize)
    if (nI != n && nI<= 1){
      testWarn="WARNING: Not enough labs were including in the consensus for DoE with LOO to be meaningful<br/>"
    }

    for (j in 1:n)
    {

      cat(j, "of", n, "\n")
      zj = linearOP(x=x[-j], u=u[-j], nu=nu[-j], weights=weights[-j], K)  ##AK change##
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
    for (j in 1:n)
    {
      cat(j, "of", n, "\n")

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

  }

  results = data.frame(Lab=lab, DoE.x=DoE.x, DoE.U95=DoE.U,
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
