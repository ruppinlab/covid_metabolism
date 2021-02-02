library(gembox)
library(stringr)
data("recon1")
nc <- 16L # number of CPUs/cores to use for parallelization

# run differential flux analysis for each dataset
# here the dataset id is passed via command line argument, such that the "dflux.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs
# alternatively, directly assign id to be the name of a particular dataset, e.g. id <- "vero"

id <- commandArgs(TRUE)
imat.res <- readRDS(paste0("imat.res.",id,".RDS"))

### run differential flux analysis:

# by each reaction
df.res <- get.diff.flux(imat.res$ctrl$result.model, imat.res$trt$result.model, nc=nc)
# total flux by each metabolite
df.met.res <- get.diff.flux.by.met(imat.res$ctrl$result.model, imat.res$trt$result.model, nc=nc)
# net transport flux between the cytosol and the extracellular space by each metabolite
df.tx.res <- get.diff.transport.flux(imat.res$ctrl$result.model, imat.res$trt$result.model, nc=nc)

df.res <- list(rxn=df.res, met=df.met.res, tx=df.tx.res)
saveRDS(df.res, file=paste0("dflux.res.",id,".RDS"))

