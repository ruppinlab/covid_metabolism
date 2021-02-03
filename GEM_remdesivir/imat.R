library(gembox)
data("recon1")
nc <- 16L # number of CPUs/cores to use for parallelization

load("data.for.gem.RData")

# run iMAT on the remdesivir treatment dataset, for each experimental group
# here the experimental group is passed via command line argument, such that the "imat.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs

x <- commandArgs(TRUE)
id <- x[1]
bm.lb <- as.numeric(x[2])

bm <- get.biomass.idx(recon1)
model <- set.rxn.bounds(recon1, bm, lbs=bm.lb, relative=TRUE)

# run iMAT

imat.res <- imat(model, exprs.int[[id]], samp.pars=list(nc=nc, n.sample=5000))

saveRDS(imat.res, file=paste0("imat.res.",id,".",bm.lb".RDS"))

