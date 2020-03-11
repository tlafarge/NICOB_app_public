######################################################################
##
## FILE        : DoEUnilateralDerSimonianLaird-AK-AP-2019Apr01.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2019-Apr-01

## INPUT : x   = Numeric vector with values measured by the labs (required)
##         u   = Numeric vector with standard uncertainties
##               associated with measured values (required)
##         nu  = Numeric vector with (positive) numbers of
##               degrees of freedom that the elements of u are
##               based on; if NULL, then the elements of u are
##               assumed to be based on infinitely many degrees
##               of freedom (optional)
##         lab = Character vector with laboratory labels (optional)
##         K   = Number of bootstrap replicates (optional)
##         LOO = TRUE/FALSE indicating whether to base the DoE on
##               differences to leave-one-out estimates of the
##               consensus value, or on differences to the
##               consensus value that has been derived from all
##               measurement results (optional)
##         DLRes = The results of the DerSimonianLaird function
##         coverageProb = Coverage probability for the interval
##               (DoE.Lwr, DoE.Upr), and for the interval DoE.x +/-
##               DoE.U (required)

## OUTPUT: List with two elements: (D) Kxn array whose jth column
## has K bootstrap replicates of the unilateral DoE for lab j;
## (DoE) data frame with values and expanded uncertainties
## for unilateral degrees of equivalence, and endpoints of coverage
## intervals for these values

######################################################################

DoEUnilateralDL = function (x, u, nu, lab, K,LOO, coverageProb, DLRes)
{
  require(metafor)
  testWarn = ""

  n = length(x)
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}

  if (LOO) {
    ## DoEs based on leave-one-out estimates
    # Check if any lab was exempt from prior Mcmc computation
    sanitize = !startsWith( lab,"-")

    #Number of labs taken into account for the consensus
    nI = sum(sanitize)
    if (nI != n && nI<= 1){
      testWarn="WARNING: Not enough labs were including in the consensus for DoE with LOO to be meaningful<br/>"
    }

    DoE.x = numeric(n)
    D = array(dim=c(K,n))
    names(DoE.x) = dimnames(D)[[2]] = lab
    for (j in 1:n)
    {
      cat (j, "of", n, "\n")
      tau2B = numeric(K)
      for (k in 1:K) {tau2B[k] = sampleFromTau2Dist(x[-j],u[-j])}

      ## AP-AK: We avoid doing the "full" bootstrap for
      ## each mu[-j], and rely instead on a modified
      ## version of the Hartung-Knapp-Sidik-Jonkman
      ## approximation for the distribution of each
      ## leave-one-out estimate of the consensus value,
      ## as described by Rover, Knapp, and Friede (2015,
      ## DOI 10.1186/s12874-015-0091-1)
      xj = x[-j]; uj = u[-j]
      nj = length(xj)
      wj0 = 1/uj^2
      xj0 = sum(wj0*xj)/sum(wj0)
      tauj2 = max(0, (sum(wj0*(xj-xj0)^2)-nj+1) /
                    (sum(wj0) - sum(wj0^2)/sum(wj0)))
      wj = 1/(uj^2+tauj2)
      muj = sum(wj*xj)/sum(wj)
      qSTAR = max(1, sum(wj*(xj-muj)^2)/(nj-1))
      muj.u = sqrt(qSTAR/sum(wj))
      muj.df = nj-1

      ## AP-AK: Addressing possibility of no
      ## more than 2 degrees of freedom
      ## Alternatives include sqrt(3) and
      ## qt((1+coverageProb)/2, df=nu[j]) /
      ## qt((1+coverageProb)/2, df=Inf)
      muj.phi = if (muj.df <= 2) {
        a = 0.4668994; b = -0.3998882
        1.5*((1 - (3/4)*(a - 4*b *
                           (muj.df^(-3/4)-1)/3))^(-4/3))
      } else {sqrt(muj.df/(muj.df-2))}

      mujB = as.vector(muj) + muj.u*rt(K, df=muj.df)/muj.phi
      DoE.x[j] = x[j] - muj
      ## AP-AK: The same sample of size K in tau2B
      ## is used for all j=1,...,n
      sigmaj = sqrt(tau2B + u[j]^2)
      if (is.null(nu[j]) || (nu[j]==Inf)) {
        ej = rnorm(K, mean=0, sd=sigmaj)
      } else {
        ## AP-AK: Addressing possibility of no
        ## more than 2 degrees of freedom (same approach as above)
        phi = if (nu[j] <= 2) {
          a = 0.4668994; b= -0.3998882
          1.5*((1 - (3/4)*(a - 4*b *
                             (nu[j]^(-3/4)-1)/3))^(-4/3))
        } else {sqrt(nu[j]/(nu[j]-2))}

        ej = sigmaj*rt(K, df=nu[j])/phi }

      ## Following Jones & Spiegelhalter (2011) --
      ## Approach 2: Identify Outliers to the Random
      ## Effects Distribution.  Note that D[,j] is used
      ## only to evaluate the uncertainty of (not to
      ## assign a value to) the difference in DoE.x
      D[,j] = x[j] + ej - mujB
    }

  } else {
    ## DoEs computed according to MRA

    z = DLRes
    mu = z$b
    indexInf = (nu == Inf)
    ## Column j of D will have a sample from the
    ## bootstrap distribution of x[j]-mu
    D = array(dim=c(K,n))
    dimnames(D)[[2]] = lab
    for (k in 1:K)
    {

      tau2B = sampleFromTau2Dist(x,u)
      xB = rnorm(n, mean=mu, sd=sqrt(tau2B + u^2))
      uB = rep(NA, n)
      if (is.null(nu)) {uB = u
      } else {
        uB[indexInf] = u[indexInf]
        uB[!indexInf] = u[!indexInf] *
          sqrt(nu[!indexInf] /
                 rchisq(sum(!indexInf), df=nu[!indexInf]))}
      muB = rma(yi=xB, sei=uB, slab=lab, method="DL",
                knha=FALSE)$b
      D[k,] = xB - as.vector(muB)
    }
    DoE.x = x - as.vector(mu)
    names(DoE.x) = lab
    ## AP 2019-Mar-30: First center each column of D at zero, then
    ## shift it to become centered at the corresponding value of DoE.x
    ## Needs to be done only for the MRA version of the DoEs
    DC = sweep(D, 2, apply(D, 2, mean))
    D = sweep(DC, 2, -DoE.x)
  }
  DoE.Lwr = apply(D, 2, function (x) {quantile(x, probs=(1-coverageProb)/2)})
  DoE.Upr = apply(D, 2, function (x) {quantile(x, probs=(1+coverageProb)/2)})
  DoE.U = rep(NA,n)
  for (j in 1:n) {
    DoE.U[j] = symmetricalBootstrapCI(D[,j], DoE.x[j], coverageProb) }

  results = data.frame(Lab=lab, DoE.x=DoE.x, DoE.U95=DoE.U,
                       DoE.Lwr=DoE.Lwr, DoE.Upr=DoE.Upr)
  return(invisible(list(D=D, DoE=results, DoEwarn=testWarn)))
}

## ######################################################################
## ##
## ## EXAMPLE
##
## source("symmetricalBootstrapCI.R")
## source("sampleFromTau2Dist.R")
##
## PCB28 = data.frame(lab=c("IRMM", "KRISS", "NARL", "NIST", "NMIJ", "NRC"),
##     x=c(34.3, 32.9, 34.53, 32.42, 31.9, 35.8),
##     u=c(1.03, 0.69, 0.83, 0.29, 0.4, 0.38),
##     nu=c(60, 4, 18, 2, 13, 60))
##
## PCB28.DL.DoEUni.MRA =
##     DoEUnilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##                     lab=PCB28$lab, K=10000, LOO=FALSE, coverageProb=0.95)
## PCB28.DL.DoEUni.MRA$DoE
## ##         Lab      DoE.x  DoE.U95   DoE.Lwr  DoE.Upr
## ## IRMM   IRMM  0.6995674 3.764803 -3.060573 4.479257
## ## KRISS KRISS -0.7004326 3.476731 -4.232422 2.678335
## ## NARL   NARL  0.9295674 3.598871 -2.706075 4.514341
## ## NIST   NIST -1.1804326 3.360617 -4.541402 2.180905
## ## NMIJ   NMIJ -1.7004326 3.289815 -4.951760 1.620506
## ## NRC     NRC  2.1995674 3.285051 -1.081067 5.502881
##
## PCB28.DL.DoEUni.LOO =
##     DoEUnilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##                     lab=PCB28$lab, K=10000, LOO=TRUE, coverageProb=0.95)
## PCB28.DL.DoEUni.LOO$DoE
## ##         Lab      DoE.x  DoE.U95   DoE.Lwr  DoE.Upr
## ## IRMM   IRMM  0.8115538 4.477546 -3.641065 5.309114
## ## KRISS KRISS -0.8426591 4.233570 -4.948352 3.617896
## ## NARL   NARL  1.0964483 4.427031 -3.259193 5.597836
## ## NIST   NIST -1.4509118 4.717553 -6.094094 3.311249
## ## NMIJ   NMIJ -2.0743890 4.272962 -6.367484 2.168603
## ## NRC     NRC  2.9009041 2.018693  0.839509 4.872593
##
## DoE.MRA = PCB28.DL.DoEUni.MRA$DoE
## DoE.LOO = PCB28.DL.DoEUni.LOO$DoE
##
## yRANGE = range(c(DoE.MRA[, "DoE.x"]-DoE.MRA[, "DoE.U95"],
##                  DoE.MRA[, "DoE.x"]+DoE.MRA[, "DoE.U95"],
##                  DoE.LOO[, "DoE.x"]-DoE.LOO[, "DoE.U95"],
##                  DoE.LOO[, "DoE.x"]+DoE.LOO[, "DoE.U95"]))
## n = nrow(DoE.MRA)
## plot(c(1,n+0.125), yRANGE, type="n", axes=FALSE,
##      xlab="", ylab="DoE / (ng/g)")
## axis(2)
## mtext(DoE.MRA[seq(1, n, 2), "Lab"], side=1, at=seq(1,n,2), line=0, cex=1)
## mtext(DoE.MRA[seq(2, n, 2), "Lab"], side=1, at=seq(2,n,2), line=1.2, cex=1)
## abline(h=0, col="Gray")
## segments(1:n, DoE.MRA[, "DoE.x"] - DoE.MRA[, "DoE.U95"],
##          1:n, DoE.MRA[, "DoE.x"] + DoE.MRA[, "DoE.U95"], lwd=2, col="Blue")
## points(1:n, DoE.MRA[, "DoE.x"], pch=18, cex=1.5, col="Red")
## segments((1:n)+0.25, DoE.LOO[, "DoE.Lwr"],
##          (1:n)+0.25, DoE.LOO[, "DoE.Upr"], lwd=2, col="DodgerBlue")
## points((1:n)+0.25, DoE.LOO[, "DoE.x"], pch=18, cex=1.5, col="Orange")
##
#####################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################
