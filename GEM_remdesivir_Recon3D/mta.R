library(gembox)
nc <- 16L # number of CPUs/cores to use for parallelization
load("../data/Recon3D.RData")
model <- recon3d

load("data.for.gem.RData")

# run rMTA on the remdesivir treatment dataset, for virus+remdesivir --> control

bm <- get.biomass.idx(model, rgx="biomass_reaction")
model <- set.rxn.bounds(model, bm, lbs=0.05, relative=TRUE)

# run rMTA

imat.res <- readRDS(file="imat.res.virus_drug.0.05.RDS")
vsrc <- rowMeans(imat.res$result.model$sample$pnts[, 4001:1e4])
rm(imat.res); gc()
mta.res <- rmta(model, vsrc, dflux.int, nc=nc, detail=FALSE)

saveRDS(mta.res, file="mta.res.RDS")

