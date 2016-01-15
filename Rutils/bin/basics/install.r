#!/usr/bin/env r
#
# A second example to install one or more packages, now with option parsing
#
# Copyright (C) 2011 - 2014  Dirk Eddelbuettel
# Copyright (C) 2014         Carl Boettiger and Dirk Eddelbuettel
# Copyright (C) 2016         Wanding Zhou
#
# Released under GPL (>= 2)

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN

doc <- "Usage: install.r [-r REPO] [-l LIBLOC] [-h] [-d DEPS] [--error] [PACKAGES ...]

-r --repos REPO     repository to install from [default: http://cran.rstudio.com]
-l --libloc LIBLOC  location in which to install [default: NA]
-d --deps DEPS      Install suggested dependencies as well [default: NA]
-e --error          Throw error and halt instead of a warning [default: FALSE]
-h --help           show this help text"

opt <- docopt(doc)

if (opt$deps == "TRUE" || opt$deps == "FALSE") {
    opt$deps <- as.logical(opt$deps)
} else if (opt$deps == "NA") {
    opt$deps <- NA
}

if (opt$libloc == "NA") opt$libloc <- NULL

if (opt$error) {
    withCallingHandlers(
        install.packages(pkgs  = opt$PACKAGES,
                         lib   = opt$libloc,
                         repos = opt$repos,
                         dependencies=opt$deps),
        warning = stop)
    
} else {
    install.packages(pkgs  = opt$PACKAGES,
                     lib   = opt$libloc,
                     repos = opt$repos,
                     dependencies=opt$deps)
}
