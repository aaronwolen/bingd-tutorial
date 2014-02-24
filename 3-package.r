library(GenomicRanges)



# GWAS object -------------------------------------------------------------

# GWAS class is a simple extension of a GRanges object with required columns
# example(.validGWAS() in GWAS-class.r)
gwas.obj

bingd:::.validGWAS(gwas.obj)

gwas.bad <- gwas.obj
mcols(gwas.bad) <- subset(mcols(gwas.bad), select = -pvalue)
bingd:::.validGWAS(gwas.bad)
rm(gwas.bad)

# as.GWAS ensures all necessary columns are present
# calculates a z-score if no zscore column is provided

# Accessors:
marker(gwas.obj)
pvalue(gwas.obj)
zscore(gwas.obj)

marker(gwas.obj)[which.min(pvalue(gwas.obj))]

# Features ----------------------------------------------------------------

# Adding multiple features to a query
f.query <- list(DNaseI = c("DNase", "k562", "k562HotSpots"),
               H3k4me3 = c("H3K4me3", "k562", "stdHotSpots"))

f.list <- hub.features(f.query, online = TRUE, genome = "hg19")

# Structure
View(f.list$DNaseI)

# Accessors 
bingd:::Cached(f.list)
bingd:::LocalPath(f.list)


# Validation
bingd:::.validFeatureList(f.list)

f.list.bad <- endoapply(f.list, subset, select = -LocalPath)
View(f.list.bad$DNaseI)

bingd:::.validFeatureList(f.list.bad)

# Caching
system(paste("open", bingd:::cache.path()))

f.list <- cache.features(f.list)


# Feature overlap ---------------------------------------------------------

?featureOverlaps

f.overlaps <- featureOverlaps(gwas.obj, f.list)

# Respects structure of FeatureList
f.overlaps
f.overlaps[[1]]
f.overlaps[[2]]


# Annotation --------------------------------------------------------------

gwas.annot <- annotate.gwas(gwas.obj, f.list)

# Features are stored as regular metadata columns (not kept in different slots)
View(mcols(gwas.annot))


# Feature groupings are maintained by adding a featureIndex slot:
slotNames(gwas.obj)
slotNames(gwas.annot)

# featureIndex provides the pointers used to retrieve features by group
f.index <- bingd:::featureIndex(gwas.annot)
f.index[[1]]
f.index[[2]]

# featureIndex accessor: fcols
# analogous to mcols, returns a flat DataFrame, ncols = # of features
fcols(gwas.annot)

# featureIndex accessor: features
# returns a list of DataFrames that respects feature groupings
# (example: consolidate() in AnnotatedGWAS-methods.r)
features(gwas.annot)


# stats functions ---------------------------------------------------------

# calc.bayes is the metafunction that calls the individual steps of analysis

# 1. calculate prior probabilities
calc.priors(gwas.obj, adjust = 2.79)

# 2. calculate conditional probabilities
calc.conditionals(gwas.annot, risk.thresh = 1e-5, adjust = 2.79)

gwas.cons <- consolidate(gwas.annot)
calc.conditionals(gwas.cons, risk.thresh = 1e-5, adjust = 2.79)

# 3. calculate posterior probabilities
# (this code should be moved to a separate function)


# automated testing -------------------------------------------------------


# furthering developent ---------------------------------------------------

# 1. using rstudio and devtools
# 2. submitting changes via github pull requests


