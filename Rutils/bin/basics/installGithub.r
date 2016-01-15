#!/usr/bin/env r
#
# A simple example to install one or more packages from GitHub
#
# Copyright (C) 2014         Carl Boettiger and Dirk Eddelbuettel
#
# Released under GPL (>= 2)

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN
suppressMessages(library(devtools))

doc <- "Usage: installGithub.r [-r REPO] [-l LIBLOC] [-h] [-d DEPS] [PACKAGES ...]

-r --repos REPO     repository to install from [default: http://cran.rstudio.com]
-l --libloc LIBLOC  location in which to install [default: /usr/local/lib/R/site-library]
-d --deps DEPS      Install suggested dependencies as well? [default: NA]
-h --help           show this help text"

opt <- docopt(doc)
if(opt$deps == "TRUE" || opt$deps == "FALSE")
    opt$deps <- as.logical(opt$deps)
if(opt$deps == "NA")
    opt$deps <- NA

options(repos = opt$repos)
install_github(repo  = opt$PACKAGES,
               paste("-l =", opt$libloc),
               dependencies = opt$deps)
