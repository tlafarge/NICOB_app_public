######################################################################
##
## FILE        : bayesLaplace.R
## AUTHOR      : Antonio Possolo
## MODIFICATION: 2018-Jan-28
##
######################################################################

## INPUT ------------------------------------------------------

##   x   = Numeric vector with n values measured by the labs
##   u   = Numeric vector with n standard uncertainties associated with
##         measured values
##   nu  = Numeric vector with n (positive or possibly +Inf) numbers of
##         degrees of freedom that the elements of u are based on; if
##         NULL, then the elements of u are assumed to be based on
##         infinitely many degrees of freedom and sigma=u
##   lab = Character vector with n laboratory labels
##   tauPriorScale
##       = Scale of the half-Cauchy prior distribution for the
##         between-lab standard deviation ("dark uncertainty") tau
##   sigmaPriorScale
##       = Scale of the half-Cauchy prior distributions for the true
##         values of the within-lab standard uncertainties sigma
##   ni  = Number of total MCMC draws for each chain
##   nb  = Number of MCMC draws to discard at the beginning of each chain
##   nt  = Thinning rate (positive integer) for each chain
##   nc  = Number of chains must be 2 or greater; if a value less than
##         2 is specified, nc is set equal to 2 (Default 2)
##   nCores
##       = Number of clusters or cores for parallel execution (Default
##         NULL means no parallelization)
##   printResults
##       = Whether to print summary statistics at the end of the
##         execution (Default TRUE)
##   progress.bar = Default "text" ("none" to suppress)
##   tempDir      = Full path to directory for temporary JAGS code

## OUTPUT -----------------------------------------------------

## Invisible list with two elements: (summary) summary
## statistics for all the parameters in the model; (mcmc) object
## of class "mcmc" as defined in R package "coda", which may
## subsequently be summarized by the "summary" function and visualized
## with the "plot" function

## MODEL & DETAILS --------------------------------------------

## The value x[j] measured by lab j=1,...,n, with associated standard
## uncertainty u[j], is represented as x[j] = xi[j] + epsilon[j].

## The lab means {xi[j]} are a sample from a Laplace distribution with
## mean mu and standard deviation tau=1/sqrt(alpha). mu is the true
## value of the consensus value.

## The lab-specific measurement errors {epsilon[j]} are independent
## Gaussian random variables all with mean 0 and standard deviations
## {sigma[j]}, where sigma[j]=u[j] when the latter is based on
## infinitely many degrees of freedom.

## When the number nu[j] of degrees of freedom that u[j] is based on
## is finite, nu[j]*(u[j]/sigma[j])^2 is modeled as a chi-squared
## random variable with nu[j] degrees of freedom.

## mu has a Gaussian prior distribution with large standard deviation
## tau has a half-Cauchy prior distribution with median tauPriorScale
## sigma[j] has a half-Cauchy prior distribution with median sigmaPriorScale

## NOTE: Precision is the reciprocal of the variance

######################################################################

bayesLaplace = function (x, u, nu=NULL,
    tauPriorScale, sigmaPriorScale,
    ni, nb, nt, nc=1,
    nCores=NULL)
{
    require(R2jags)
    n = length(x)

    ## In JAGS (version 4.3.0), the Laplace distribution with mean mu
    ## and precision alpha is defined by y ~ ddexp(mu, alpha). Its
    ## probability density function is alpha*exp(-alpha*abs(y-mu))/2,
    ## and its variance is 2/(alpha ^2)

   if (is.null(nu)) {
       ## When nu is NULL we take sigma=u, as if all
       ## elements of u were based on infinitely many
       ## degrees of freedom, and the sigma are not estimated
       cat("model {
                   for(j in 1:n)
                     {
                       xi[j] ~ ddexp(mu, alpha)
                       sigmaPrecision[j] <- pow(u[j], -2)
                       x[j] ~ dnorm(xi[j], sigmaPrecision[j])
                     }
                   mu ~ dnorm(0, 1E-10)
                   tauPriorPrecision <- pow(tauPriorScale, -2)
                   tau ~ dt(0, tauPriorPrecision, 1)T(0,)
                   alpha <- sqrt(2)/tau
                  }", fill=TRUE, file="LabComp-Model.txt")
       LabComp.Data = list(n=n, x=x, u=u, tauPriorScale=tauPriorScale)

       ## Inits function
       LabComp.inits = function() {
           list(mu = rnorm(1, mean=median(x), sd=mad(x)),
                tau = mad(x),
                xi = rnorm(n, mean=median(x), sd=u) ) }

       ## Parameters to estimate
       LabComp.params = c("mu", "tau", "xi")
   } else {
       ## When nu is not NULL and there are some elements of
       ## nu that are Inf, then we replace Inf by 1000
       nu[nu == Inf] = 1000
       cat("model {
                   for(j in 1:n)
                     {
                       xi[j] ~ ddexp(mu, alpha)
                       sigma[j] ~ dt(0, sigmaPriorPrecision, 1)T(0,)
                       sigmaPrecision[j] <- pow(sigma[j], -2)
                       x[j] ~ dnorm(xi[j], sigmaPrecision[j])
                       varmod[j]~ dgamma(nu[j]/2, 1/(2*sigma[j]*sigma[j]))
                     }
                   sigmaPriorPrecision <- pow(sigmaPriorScale, -2)
                   mu ~ dnorm(0, 1E-10)
                   tauPriorPrecision <- pow(tauPriorScale, -2)
                   tau ~ dt(0, tauPriorPrecision, 1)T(0,)
                   alpha <- sqrt(2)/tau
                  }", fill=TRUE, file="LabComp-Model.txt")
       LabComp.Data = list(n=n, x=x, nu=nu, varmod=nu*u^2,
                           tauPriorScale=tauPriorScale,
                           sigmaPriorScale=sigmaPriorScale )

        ## Inits function
        LabComp.inits = function() {
            list(mu = rnorm(1, mean=median(x), sd=mad(x)),
                 tau = mad(x),
                 xi = rnorm(n, mean=median(x), sd=u),
                 sigma = u) }

        ## Parameters to estimate
        LabComp.params = c("mu", "tau", "xi", "sigma")
    }

    ## Start Gibbs sampling

    withProgress(message = 'Running Bayes', style="old",value = 0, {
      LabComp.output =
        if (is.null(nCores)) {
            jags(LabComp.Data, inits=LabComp.inits, parameters.to.save=LabComp.params,
                 model.file="LabComp-Model.txt", n.thin=nt, n.chains=nc,
                 n.burnin=nb, n.iter=ni)
        } else { #leaving this in as a potential future direction, NICOB currently does not allow for this.
            jags.parallel(LabComp.Data, inits=LabComp.inits,
                          parameters.to.save=LabComp.params,
                          model.file="LabComp-Model.txt",
                          nCores=min(c(nCores, nc)),
                          n.thin=nt, n.chains=nc,
                          n.burnin=nb, n.iter=ni) }
    })
    bayesres = as.mcmc(LabComp.output)
    geweke = geweke.diag(bayesres)
    zScores = geweke[[1]][[1]][-1] # confirm that this index is correct
    # for (jc in 2:nc) {zScores = c(zScores, geweke[[jc]][[1]][-1])} #changed nc to 1
    geweke.pvalues = p.adjust(2*pnorm(-abs(zScores)), method="BH")
    if (min(geweke.pvalues) < 0.05) {
        cat("WARNING: MCMC may not have reached equilibrium\n")
        cat("WARNING: Results are unreliable\n")
        cat(paste("SUGGESTION: Re-run with ",
                  "ni = ", 10*ni, ", nb = ", 2*ni,
                  ", and nt = ", ceiling(0.001*ni), "\n", sep=""))

        testWarn=paste("WARNING: MCMC may not have reached equilibrium<br/>",
                       "WARNING: Results are unreliable<br/>",
                       paste("SUGGESTION: Re-run with ",
                             "'Total number of iterations' = ", 10*ni, ", 'Length of burn in' = ", 2*nb,
                             ", and 'Thinning rate' = ", ceiling(0.001*ni), "<br/>", sep=""),sep=" ")
    }else{
      testWarn=""
    }

    if (file.exists("LabComp-Model.txt")) {file.remove("LabComp-Model.txt")}

    return(list(mcmcout=bayesres,warn=testWarn))
}
