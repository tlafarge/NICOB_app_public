######################################################################
##
## FILE        : DoEBilateralBayes.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2016-Jul-05
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
##         LOO = TRUE/FALSE indicating whether to base the DoE on
##               differences to leave-one-out estimates of the
##               consensus value, or on differences to the
##               consensus value that has been derived from all
##               measurement results (optional)
##         DoEUnilateral = Output of DoEUnilateralBayes (optional)

## OUTPUT: List with four elements: (B.x) nxn matrix with the
## estimates of the between-lab differences; (B.U95), nxn matrix
## with the corresponding expanded uncertainties for 95 %
## coverage; (B.Lwr95) nxn matrix with the 2.5th percentiles of
## the pairwise differences; (B.Upr95) nxn matrix with the 97.5th
## percentiles of the pairwise differences.

######################################################################

DoEBilateralBayes = function (x, u, nu, lab, #LOO,
                              DoEUnilateral,coverageProb)
{
  n = length(x)
  if (is.null(lab)) {lab=paste("L", 1:n, sep="")}
  sanitize = !startsWith( lab,"-")

  #Number of labs taken into account for the consensus
  nI = sum(sanitize)
  if (nI<= 2) {
    testWarn="WARNING: Not enough labs were including in the consensus for DoE to be available with this method<br/>"
  }

  # DoEUni = if (is.null(DoEUnilateral)) {
  #   DoEUnilateralBayes(x, u, nu, lab, LOO=LOO) } else {
  #     DoEUnilateral}
  DoEUni = DoEUnilateral

  DoE.Bilateral.x = DoE.Bilateral.U = diag(rep(0,n))
  DoE.Bilateral.Lwr = DoE.Bilateral.Upr = diag(rep(0,n))
  dimnames(DoE.Bilateral.x) = list(lab, lab)
  dimnames(DoE.Bilateral.U) = list(lab, lab)
  dimnames(DoE.Bilateral.Lwr) = list(lab, lab)
  dimnames(DoE.Bilateral.Upr) = list(lab, lab)

  #DoE not possible if we have less than 3 labs in the consensus
  if (nI == n || nI >= 3) {
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
  }
  return(list(B.x=DoE.Bilateral.x,
              B.U=DoE.Bilateral.U,
              B.Lwr=DoE.Bilateral.Lwr,
              B.Upr=DoE.Bilateral.Upr))
}
