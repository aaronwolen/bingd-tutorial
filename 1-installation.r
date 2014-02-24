# bingd tutorial

install.packages(c("devtools", "plyr"))
library(devtools)

bioc.deps <- c("GenomicRanges", "AnnotationHub")

source("http://www.bioconductor.org/biocLite.R")
biocLite(bioc.deps, suppressUpdates = TRUE)


gh.username <- "username"
gh.password <- "password"

install_github("bingd", "aaronwolen", 
               auth_user = gh.username, 
               password = gh.password)

# data url: http://bit.ly/bingd-files
