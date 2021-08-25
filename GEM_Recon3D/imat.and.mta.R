library(gembox)
nc <- 16L # number of CPUs/cores to use for parallelization
load("../data/Recon3D.RData")
model <- recon3d

dat <- readRDS("data.for.gem.RDS")

# run iMAT and rMTA on each dataset
# here the dataset id is passed via command line argument, such that the "imat.and.mta.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs

x <- commandArgs(TRUE)
id <- x[1] # dataset ID
bm.lb <- as.numeric(x[2]) # minimal biomass requirement

exprs.int <- dat$exprs.int[[id]]

if (bm.lb!=0) {
  bm <- get.biomass.idx(model, rgx="biomass_reaction")
  model <- set.rxn.bounds(model, bm, lbs=bm.lb, relative=TRUE)
}

# run iMAT

ctrl <- imat(model, exprs.int$ctrl, samp.pars=list(nc=nc, n.sample=1e4)) # control samples
trt <- imat(model, exprs.int$trt, samp.pars=list(nc=nc, n.sample=1e4)) # virus-infected samples
imat.res <- list(ctrl=ctrl, trt=trt)

saveRDS(imat.res, file=paste0("imat.res.",id,".RDS"))

# run rMTA

vtrt <- rowMeans(imat.res$trt$result.model$sample$pnts[, 4001:1e4])
mta.res <- rmta(model, vtrt, dflux.int, nc=nc, detail=FALSE)

saveRDS(mta.res, file=paste0("mta.res.",id,".RDS"))