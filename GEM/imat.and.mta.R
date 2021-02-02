library(gembox)
data("recon1")
data("media")
nc <- 16L # number of CPUs/cores to use for parallelization
model <- recon1

dat <- readRDS("data.for.gem.RDS")

# run iMAT and rMTA on each dataset
# here the dataset id is passed via command line argument, such that the "imat.and.mta.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs

x <- commandArgs(TRUE)
id <- x[1] # dataset ID
media <- x[2] # media
bm.lb <- as.numeric(x[3]) # minimal biomass requirement

exprs.int <- dat$exprs.int[[id]]
dflux.int <- dat$dflux.int[[id]]

if (media=="1") model <- set.medium(model, media$dmem.lo.glc, set.all=TRUE)
if (bm.lb!=0) {
  bm <- get.biomass.idx(model)
  model <- set.rxn.bounds(model, bm, lbs=bm.lb, relative=TRUE)
}

# run iMAT

ctrl <- imat(model, exprs.int$ctrl, samp.pars=list(nc=nc, n.sample=5000)) # control samples
trt <- imat(model, exprs.int$trt, samp.pars=list(nc=nc, n.sample=5000)) # virus-infected samples
imat.res <- list(ctrl=ctrl, trt=trt)

saveRDS(imat.res, file=paste0("imat.res.",id,".RDS"))

# run rMTA

vtrt <- rowMeans(imat.res$trt$result.model$sample$pnts[, 2001:5000])
mta.res <- rmetal(model, vtrt, dflux.int, nc=nc, detail=FALSE)

saveRDS(mta.res, file=paste0("mta.res.",id,".RDS"))

