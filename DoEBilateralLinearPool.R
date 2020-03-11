######################################################################
##
## FILE        : DoEBilateralLinearPool.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Jun-05
##
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
##
## OUTPUT: List with four elements: (B.x) nxn matrix with the
## estimates of the between-lab differences; (B.U95), nxn matrix
## with the corresponding expanded uncertainties for 95 %
## coverage; (B.Lwr95) nxn matrix with the 2.5th percentiles of
## the pairwise differences; (B.Upr95) nxn matrix with the 97.5th
## percentiles of the pairwise differences.

######################################################################

DoEBilateralPool = function (x, u, nu, lab, weights, K,coverageProb) ##AK change##
{
  source("linearOP.R") ##AK change##
  
  n = length(x)
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}
  
  # z28.POOL = linearPool(x=x, u=u, lab=lab, nu=nu) ##AK change## this isn't used later
  DoE.Bilateral.POOL.x = DoE.Bilateral.POOL.U = diag(rep(0,n))
  DoE.Bilateral.POOL.Lwr = DoE.Bilateral.POOL.Upr = diag(rep(0,n))
  dimnames(DoE.Bilateral.POOL.x) = list(lab, lab)
  dimnames(DoE.Bilateral.POOL.U) = list(lab, lab)
  dimnames(DoE.Bilateral.POOL.Lwr) = list(lab, lab)
  dimnames(DoE.Bilateral.POOL.Upr) = list(lab, lab)
  for (j1 in 1:(n-1))
  {
    cat(j1, "of", n-1, "\n")
    zj1 = linearOP(x=x[-j1], u=u[-j1], nu=nu[-j1],weights=weights[-j1],K) ##AK change##
    ej1 = sample(zj1-mean(zj1), size=K, replace=TRUE)
    Dj1 = x[j1] + ej1 - zj1
    for (j2 in (j1+1):n)
    {
      zj2 = linearOP(x=x[-j2], u=u[-j2], nu=nu[-j2],weights=weights[-j2],K) ##AK change##
      ej2 = sample(zj2-mean(zj2), size=K, replace=TRUE)
      Dj2 = x[j2] + ej2 - zj2
      DoE.Bilateral.POOL.x[j1,j2] = x[j1]-x[j2]
      DoE.Bilateral.POOL.x[j2,j1] = DoE.Bilateral.POOL.x[j1,j2] 
      # DoE.Bilateral.POOL.x[j2,j1] = x[j2]-x[j1]
      DoE.Bilateral.POOL.U[j1,j2] =
        symmetricalBootstrapCI(Dj1-Dj2, x[j1]-x[j2], coverageProb)
      DoE.Bilateral.POOL.U[j2,j1] =
        symmetricalBootstrapCI(Dj2-Dj1, x[j2]-x[j1], coverageProb)
      DoE.Bilateral.POOL.Lwr[j1,j2] =
        quantile(Dj1-Dj2, probs=(1-coverageProb)/2)#0.025)
      DoE.Bilateral.POOL.Lwr[j2,j1] =
        quantile(Dj2-Dj1, probs=(1-coverageProb)/2)#0.025)
      DoE.Bilateral.POOL.Upr[j1,j2] =
        quantile(Dj1-Dj2, probs=(1+coverageProb)/2)#0.975)
      DoE.Bilateral.POOL.Upr[j2,j1] =
        quantile(Dj2-Dj1, probs=(1+coverageProb)/2)#0.975)
    }
  }
  
  return(list(B.x=DoE.Bilateral.POOL.x,
              B.U=DoE.Bilateral.POOL.U,
              B.Lwr=DoE.Bilateral.POOL.Lwr,
              B.Upr=DoE.Bilateral.POOL.Upr))
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
## DoE = DoEBilateralPool(x=PCB28$x, u=PCB28$u, nu=PCB28$nu, lab=PCB28$lab)
##
## round(DoE$B.x, 2)
## ##        IRMM KRISS  NARL  NIST NMIJ   NRC
## ## IRMM   0.00  1.40 -0.23  1.88 2.40 -1.50
## ## KRISS -1.40  0.00 -1.63  0.48 1.00 -2.90
## ## NARL   0.23  1.63  0.00  2.11 2.63 -1.27
## ## NIST  -1.88 -0.48 -2.11  0.00 0.52 -3.38
## ## NMIJ  -2.40 -1.00 -2.63 -0.52 0.00 -3.90
## ## NRC    1.50  2.90  1.27  3.38 3.90  0.00
##
## round(DoE$B.U95, 2)
## ##       IRMM KRISS NARL NIST NMIJ  NRC
## ## IRMM  0.00  6.16 6.10 6.05 5.87 5.55
## ## KRISS 6.16  0.00 6.17 6.04 5.83 5.66
## ## NARL  6.10  6.17 0.00 6.07 5.88 5.54
## ## NIST  6.05  6.04 6.07 0.00 5.69 5.59
## ## NMIJ  5.87  5.83 5.88 5.69 0.00 5.41
## ## NRC   5.55  5.66 5.54 5.59 5.41 0.00
##
## DoE.Bilateral.LinearPool.x = DoE$B.x
## DoE.Bilateral.LinearPool.U95 = DoE$B.U95
## xx = ((((DoE.Bilateral.LinearPool.x - DoE.Bilateral.LinearPool.U95) *
##             (DoE.Bilateral.LinearPool.x +
##                  DoE.Bilateral.LinearPool.U95)) < 0))
## xx[xx == 1] = 0.9
## n = dim(DoE.Bilateral.LinearPool.x)[[1]]
## diag(xx) = rep(1, n)
## image(1:n, 1:n, xx, xlim=c(-1, n+1),
##       breaks=c(-0.5, 0.5, 0.95, 1.5), col=heat.colors(3),
##       axes=FALSE, xlab="", ylab="", asp=1)
## polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
##         border="Blue", lwd=2)
## text(rep(0, n), 1:n, PCB28$lab, adj=1, col="Blue", font=2)
## mtext(PCB28$lab, side=1, at=1:n, line=-0.75, adj=0.5, col="Blue", font=2)
## mtext("PCB 28 Bilateral DoE: Linear Pool", side=3, at=0, line=1.5,
##       adj=0, col="Blue", font=2)

######################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################
