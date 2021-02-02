library(gembox)
data("recon1")
nc <- 16L # number of CPUs/cores to use for parallelization

load("data.for.gem.RData")

# run rMTA on the remdesivir treatment dataset, between selected pairs of groups
# here the id (representing pairs of groups) is passed via command line argument, such that the "mta.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs
id <- commandArgs(TRUE)

bm <- get.biomass.idx(recon1)
model <- set.rxn.bounds(recon1, bm, lbs=0.05, relative=TRUE)

if (id %in% c("v2c","v2vd")) {
  grp <- "virus"
} else if (id %in% c("vd2c","vd2d")) {
  grp <- "virus_drug"
} else if (id=="d2c") grp <- "drug"

# run rMTA

imat.res <- readRDS(file=paste0("imat.res.",grp,".RDS"))
vsrc <- rowMeans(imat.res$result.model$sample$pnts[, 2001:5000])

mta.res <- rmetal(model, vsrc, dflux.int[[id]], nc=nc, detail=FALSE)

saveRDS(mta.res, file=paste0("mta.res.",id,".RDS"))

