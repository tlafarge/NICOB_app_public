######################################################################
##
## FILE        : DoEBilateralDerSimonianLaird.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Jun-26

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
##               consensus value or on differences to the
##               consensus value derived from all measurement
##               results (optional)
##         DoEUnilateral = Output of DoEUnilateralDL (optional)

## OUTPUT: List with four elements: (B.x) nxn matrix with the
## estimates of the between-lab differences; (B.U95), nxn matrix
## with the corresponding expanded uncertainties for 95 %
## coverage; (B.Lwr95) nxn matrix with the 2.5th percentiles of
## the pairwise differences; (B.Upr95) nxn matrix with the 97.5th
## percentiles of the pairwise differences.

######################################################################

DoEBilateralDL = function (x, u, nu, lab, K,
                           # LOO,
                           coverageProb,DoEUnilateral)
{

  lab = dimnames(DoEUnilateral$DoE)[[1]]
  DoEUni = DoEUnilateral#}

  n = length(x)
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}

  DoE.Bilateral.x = DoE.Bilateral.U = diag(rep(0,n))
  DoE.Bilateral.Lwr = DoE.Bilateral.Upr = diag(rep(0,n))
  dimnames(DoE.Bilateral.x) = list(lab, lab)
  dimnames(DoE.Bilateral.U) = list(lab, lab)
  dimnames(DoE.Bilateral.Lwr) = list(lab, lab)
  dimnames(DoE.Bilateral.Upr) = list(lab, lab)
  for (j1 in 1:(n-1))
  {
    for (j2 in (j1+1):n)
    {
      Bj1j2 = DoEUni$D[,j1] - DoEUni$D[,j2]
      DoE.Bilateral.x[j1,j2] =
        DoEUni$DoE[j1,"DoE.x"] - DoEUni$DoE[j2,"DoE.x"]
      DoE.Bilateral.x[j2,j1] = DoE.Bilateral.x[j1,j2]
      DoE.Bilateral.Lwr[j1,j2] = quantile(Bj1j2, probs=(1-coverageProb)/2)#0.025)
      DoE.Bilateral.Lwr[j2,j1] = DoE.Bilateral.Lwr[j1,j2]
      DoE.Bilateral.Upr[j1,j2] = quantile(Bj1j2, probs=(1+coverageProb)/2)#0.975)
      DoE.Bilateral.Upr[j2,j1] = DoE.Bilateral.Upr[j1,j2]
      DoE.Bilateral.U[j1,j2] =
        symmetricalBootstrapCI(Bj1j2, mean(Bj1j2), coverageProb)
      DoE.Bilateral.U[j2,j1] = DoE.Bilateral.U[j1,j2]
    }
  }

  return(list(B.x=DoE.Bilateral.x,
              B.U=DoE.Bilateral.U,
              B.Lwr=DoE.Bilateral.Lwr,
              B.Upr=DoE.Bilateral.Upr))
}

## ######################################################################
## ##
## ## EXAMPLE
##
## source("~/NIST/TALES/Analysis/symmetricalBootstrapCI.R")
## source("~/NIST/TALES/Analysis/DoEUnilateralDerSimonianLaird.R")
## source("~/NIST/TALES/Analysis/sampleFromTau2Dist.R")
##
## PCB28 = data.frame(lab=c("IRMM", "KRISS", "NARL", "NIST", "NMIJ", "NRC"),
##     x=c(34.3, 32.9, 34.53, 32.42, 31.9, 35.8),
##     u=c(1.03, 0.69, 0.83, 0.29, 0.4, 0.38),
##     nu=c(60, 4, 18, 2, 13, 60))
##
## ## ------------------------------------------------------------
## ## CONVENTIONAL (LOO = FALSE)
##
## ## AP-AK: I'm worried that with LOO=FALSE (MRA) we do not find any
## ## significant bilateral differences involving NRC
##
## DoEUni.MRA = DoEUnilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##     lab=PCB28$lab, LOO=FALSE, K=10000)
## DoEBi.MRA = DoEBilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##     lab=PCB28$lab, LOO=FALSE, DoEUnilateral=DoEUni.MRA)
##
## round(DoEBi.MRA$B.x, 2)
## ##        IRMM KRISS  NARL  NIST NMIJ   NRC
## ## IRMM   0.00  1.40 -0.23  1.88 2.40 -1.50
## ## KRISS -1.40  0.00 -1.63  0.48 1.00 -2.90
## ## NARL   0.23  1.63  0.00  2.11 2.63 -1.27
## ## NIST  -1.88 -0.48 -2.11  0.00 0.52 -3.38
## ## NMIJ  -2.40 -1.00 -2.63 -0.52 0.00 -3.90
## ## NRC    1.50  2.90  1.27  3.38 3.90  0.00
##
## round(DoEBi.MRA$B.U95, 2)
## ##       IRMM KRISS NARL NIST NMIJ  NRC
## ## IRMM  0.00  5.62 5.48 5.50 5.47 5.38
## ## KRISS 5.62  0.00 5.37 5.21 5.34 5.38
## ## NARL  5.48  5.37 0.00 5.26 5.21 5.29
## ## NIST  5.50  5.21 5.26 0.00 5.16 5.15
## ## NMIJ  5.47  5.34 5.21 5.16 0.00 5.12
## ## NRC   5.38  5.38 5.29 5.15 5.12 0.00
##
## DoE.Bilateral.DL.x = DoEBi.MRA$B.x
## DoE.Bilateral.DL.Lwr95 = DoEBi.MRA$B.Lwr95
## DoE.Bilateral.DL.Upr95 = DoEBi.MRA$B.Upr95
##
## n = nrow(PCB28)
## labNAMEs = dimnames(DoE.Bilateral.DL.x)[[1]]
## xx = (DoE.Bilateral.DL.Lwr95 * DoE.Bilateral.DL.Upr95 > 0)
## ij = which(xx, arr.ind=TRUE, useNames = TRUE)
## xx[xx == 1] = 0.9
## xx[xx == 0] = 0.1
## diag(xx) = rep(0.5, n)
## image(1:n, 1:n, xx, xlim=c(-1.5, n+1),
##       col=heat.colors(3),
##       axes=FALSE, xlab="", ylab="", asp=1)
## K = nrow(ij)
## if (K > 0) {for (k in 1:K) {points(ij[k,1], ij[k,2], pch=8, font=3, col="Blue")}}
## polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
##         border="Blue", lwd=2)
## text(rep(0, n), 1:n, labNAMEs, adj=1, col="Blue", font=2)
## mtext(labNAMEs, side=1, at=1:n, line=0.75, adj=0.5, col="Blue", font=2)
##
## ## ------------------------------------------------------------
## ## LEAVE-ONE-OUT (LOO = TRUE)
##
## DoEUni.LOO = DoEUnilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##     lab=PCB28$lab, LOO=TRUE, K=10000)
## DoEBi.LOO = DoEBilateralDL(x=PCB28$x, u=PCB28$u, nu=PCB28$nu,
##     lab=PCB28$lab, LOO=TRUE, DoEUnilateral=DoEUni.LOO)
##
## round(DoEBi.LOO$B.x, 2)
## ##        IRMM KRISS  NARL  NIST NMIJ   NRC
## ## IRMM   0.00  1.65 -0.28  2.26 2.89 -2.09
## ## KRISS -1.65  0.00 -1.94  0.61 1.23 -3.74
## ## NARL   0.28  1.94  0.00  2.55 3.17 -1.80
## ## NIST  -2.26 -0.61 -2.55  0.00 0.62 -4.35
## ## NMIJ  -2.89 -1.23 -3.17 -0.62 0.00 -4.98
## ## NRC    2.09  3.74  1.80  4.35 4.98  0.00
##
## round(DoEBi.LOO$B.U95, 2)
## ##       IRMM KRISS NARL NIST NMIJ  NRC
## ## IRMM  0.00  5.75 5.89 8.28 5.64 5.69
## ## KRISS 5.75  0.00 5.68 8.05 5.54 5.42
## ## NARL  5.89  5.68 0.00 8.38 5.57 5.51
## ## NIST  8.28  8.05 8.38 0.00 8.14 8.15
## ## NMIJ  5.64  5.54 5.57 8.14 0.00 5.41
## ## NRC   5.69  5.42 5.51 8.15 5.41 0.00
##
## DoE.Bilateral.DL.x = DoEBi.LOO$B.x
## DoE.Bilateral.DL.Lwr95 = DoEBi.LOO$B.Lwr95
## DoE.Bilateral.DL.Upr95 = DoEBi.LOO$B.Upr95
##
## n = nrow(PCB28)
## labNAMEs = dimnames(DoE.Bilateral.DL.x)[[1]]
## xx = (DoE.Bilateral.DL.Lwr95 * DoE.Bilateral.DL.Upr95 > 0)
## ij = which(xx, arr.ind=TRUE, useNames = TRUE)
## xx[xx == 1] = 0.9
## xx[xx == 0] = 0.1
## diag(xx) = rep(0.5, n)
## image(1:n, 1:n, xx, xlim=c(-1.5, n+1),
##       col=heat.colors(3),
##       axes=FALSE, xlab="", ylab="", asp=1)
## K = nrow(ij)
## if (K > 0) {for (k in 1:K) {points(ij[k,1], ij[k,2], pch=8, font=3, col="Blue")}}
## polygon(c(0.5,n+0.5,n+0.5,0.5), c(0.5,0.5,n+0.5,n+0.5),
##         border="Blue", lwd=2)
## text(rep(0, n), 1:n, labNAMEs, adj=1, col="Blue", font=2)
## mtext(labNAMEs, side=1, at=1:n, line=0.75, adj=0.5, col="Blue", font=2)

######################################################################
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
######################################################################
