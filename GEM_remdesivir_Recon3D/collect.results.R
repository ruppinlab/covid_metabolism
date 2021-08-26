library(stringr)

# collect the differential flux analysis results for the remdesivir treatment data
# save all results into a single data file for the convenience of downstream analysis

fns <- dir(pattern="^dflux\\.res")
names(fns) <- str_sub(fns, 11, -5)
tmp <- lapply(fns, readRDS)
df.res <- lapply(tmp, function(x) x$rxn)
df.met.res <- lapply(tmp, function(x) x$met)
df.tx.res <- lapply(tmp, function(x) x$tx)
save(df.res, df.met.res, df.tx.res, file="dflux.res.RData")

