# Introduction

NICOB.app is an R package created to make it easier to run a local
version of the NICOB app. This manual explain how to install it and all
its dependencies on a Ubuntu or any other debian based system.

# R

At the time of writing this article, the latest stable version of R is
version 3.6. We recommend running NICOB with this version.

The R packages from the Ubuntu repositories are often outdated so we
will install R by adding the repository maintained by CRAN.

To install the latest stable version of R on Ubuntu, follow these steps:
Install the packages necessary to add a new repository over HTTPS:

    sudo apt install apt-transport-https software-properties-common

Enable the CRAN repository and add the CRAN GPG key to your system using
the following commands:

    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

Now that the apt repository is added, update the packages list and
install R by typing:

    sudo apt update
    sudo apt install r-base

# System dependencies

To be able to run NICOB.app you will need to install R2jags and shiny,
two R packages that require the following system dependencies g++ and
jags. You can install them both with the commands

    sudo apt install g++
    sudo apt install jags

# R packages

To be able tun run NICOB you will need several R packages.

Launch R with:

```
R
```

And install the required packages from R with:

    install.packages("shiny")
    install.packages("metafor")
    install.packages("R2jags")

# NICOB.app

As NICOB.app is not yet released on CRAN you will need to install it
from source you can do so with the following commands:

Get the last version of the package using ```wget```, or you can use your browser to download it.

    wget https://github.com/tlafarge/NICOB_app_public/releases/download/1.2/NICOB.app_1.2.tar.gz

You can then install it with:

    R CMD INSTALL NICOB.app_1.2.tar.gz

You can now run the app with the following commands, first run R:

```
R
```

then Load the package and call the function to run the app. It will run
NICOB on a the specified port.

    require("NICOB.app")
    run_NICOB_app()

It should give you the address of the app, if you did not specify a port:

    Listening on http://127.0.0.1:3999

If it does not open your browser automatically you can navigate to this
page manually to find the app. You can also use a link like this to load
an example.

    http://127.0.0.1:3999?example=carotid
    http://127.0.0.1:3999?example=pcb
    http://127.0.0.1:3999?example=triple
