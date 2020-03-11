# NIST Consensus Builder
![NIST](https://nccoe.nist.gov/sites/all/themes/custom/nccoe2x/asset/img/NIST_logo.svg)

 The NIST Consensus Machine https://consensus.nist.gov  or NICOB
 serves to combine measurement results obtained by different laboratories or by application of different measurement methods, into a consensus estimate of the value of a scalar measurand. The NICOB qualifies the consensus estimate with an evaluation of measurement uncertainty that captures not only the stated uncertainties associated with the individual measured values, but also any additional component of uncertainty that manifests itself only when these measured values are inter-compared.


## Manual
[PDF available here](https://consensus.nist.gov/NISTConsensusBuilder-UserManual.pdf)

## Requirements
* R version 3.6+
* Jags
* R packages shiny, metafor and R2jags

## Changelog
### version 1.2
 - Added the possibility to remove labs from the consensus computation but not the degrees of equivalences by adding a "-" in front of their label

### version 1.1
  - Improve the configuration files save and load feature
  - Added the possibility to specify the wamted number of significant digits
  - Nows display the dark Uncertainty

### version 1.0
  - First public release of NICOB
