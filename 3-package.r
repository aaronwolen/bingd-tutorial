library(GenomicRanges)


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



