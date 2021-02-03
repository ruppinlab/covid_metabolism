library(gembox)
data("recon1")
nc <- 16L # number of CPUs/cores to use for parallelization

load("data.for.gem.RData")

# run rMTA on the remdesivir treatment dataset, for virus+remdesivir --> control

bm <- get.biomass.idx(recon1)
model <- set.rxn.bounds(recon1, bm, lbs=0.05, relative=TRUE)

# run rMTA

imat.res <- readRDS(file="imat.res.virus_drug.0.05.RDS")
vsrc <- rowMeans(imat.res$result.model$sample$pnts[, 2001:5000])

mta.res <- rmetal(model, vsrc, dflux.int$vd2c, nc=nc, detail=FALSE)

saveRDS(mta.res, file="mta.res.RDS")

