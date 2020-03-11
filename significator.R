######################################################################

## significator.R
##
## AUTHOR:            Thomas Lafarge & Antonio Possolo
## MODIFICATION DATE: 2017 Nov 15

######################################################################

significator = function (y, nsubsamples=1000, stat=mean, k=1,
                         digits=FALSE, tol=sqrt(.Machine$double.eps), ...)
{
    ny = length(y)
    ys = stat(y, ...)
    if (abs(ys-mean(y)) < tol) {sigma = sd(y)/sqrt(ny)
    } else {
        nz = nsubsamples
        ## MZ is the size of each sub-sample
        mz = trunc(ny/nz)
        if (mz < 25) {
            cat(paste0("\r\n",
                       "## ERROR: Number of significant digits ",
                       "cannot be computed\r\n",
                       "##        unless the size of each subsample\r\n",
                       "##        is 25 at least\r\n"))
            stop()
        }
        z = array(y, dim=c(mz,nz))
        zs = apply(z, 2, stat, ...)
        ## SIGMA is the sub-sampling estimate of the standard error of STAT(Y)
        sigma = mad(zs)/sqrt(ny/nz)
    }
    stat.signif = floor(log10(abs(ys/(k*sigma))))
    if (digits) {return(stat.signif)
    } else {
        ys.signif = signif(ys, stat.signif)
        names(ys.signif) = NULL
        return(ys.signif)
    }
}
