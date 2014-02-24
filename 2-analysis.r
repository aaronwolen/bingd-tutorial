library(bingd)

library(GenomicRanges)
library(bingdviz)

# Create GWAS object ------------------------------------------------------

gwas.df <- read.csv("data/pgc-scz-hg19.csv.gz")


gwas.obj <- as.GWAS(gwas.df, genome = "hg19", 
                    marker="snp", chr="chr", bp="bp", 
                    pvalue="pval", or="or")


# Identify features -------------------------------------------------------

# FeatureList for Conservation from local file
cons.query <- list(Conservation = "phastcon")
cons.list <- local.features(cons.query, path = "data")

# FeatureList for DNaseI data from AnnotationHub
dnase.query <- list(DNaseI = c("DNase", "frontal|neurosphere"))
dnase.list <- hub.features(dnase.query, online = TRUE, genome = "hg19")

# Download uncached features
dnase.list <- cache.features(dnase.list)


# Combine FeatureLists
feature.list <- c(dnase.list, cons.list)


# Feature enrichment ------------------------------------------------------

gwas.obj$log.pvalue <- -log10(pvalue(gwas.obj))
thresh.levels <- seq(1, 11, 1)

enrich.df <- calc.enrich(gwas.obj, feature.list,
                         stat = gwas.obj$log.pvalue, 
                         thresh.levels = thresh.levels)

View(enrich.df)

enrich.plot(enrich.df, count.thres = 3)


# GWAS annotation ---------------------------------------------------------

?annotate.gwas
gwas.annot <- annotate.gwas(gwas.obj, feature.list)

summary(gwas.annot)
summary(gwas.annot[pvalue(gwas.annot) < 1e-6])

enrich.df <- calc.enrich(gwas.annot, 
                         stat = gwas.annot$log.pvalue, 
                         thresh.levels = thresh.levels)
enrich.plot(enrich.df)


# Bayes analysis ----------------------------------------------------------

?consolidate
gwas.annot <- consolidate(gwas.annot)

summary(gwas.annot)

# Try with and without verbose = TRUE
gwas.probs <- calc.bayes(gwas.annot, risk.thresh = 1e-5, adjust = 2.79)


# Regions of interest -----------------------------------------------------
locus <- GRanges("chr7", IRanges(104e6, 106e6))

region.plot(gwas.probs, range = locus)

region.plot(gwas.probs, range = locus, y = "log.pvalue", y.axis = "-log10 p-value")



