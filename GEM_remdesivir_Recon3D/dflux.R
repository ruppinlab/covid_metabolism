library(gembox)
nc <- 16L # number of CPUs/cores to use for parallelization
load("../data/Recon3D.RData")
model <- recon3d

# run differential flux analysis for the remdesivir treatment dataset between pairs of groups
# here the id (representing pairs of groups) is passed via command line argument, such that the "dflux.swarm" in this folder can be used on a HPC server to generate an array/swarm of jobs
id <- commandArgs(TRUE)
g1 <- switch(id, c2v=, c2vd=, c2d="ctrl", v2vd="virus", d2vd="drug")
g2 <- switch(id, c2vd=, v2vd=, d2vd="virus_drug", c2v="virus", c2d="drug")

imat.res1 <- readRDS(paste0("imat.res.",g1,".0.RDS"))
imat.res2 <- readRDS(paste0("imat.res.",g2,".0.RDS"))

### run differential flux analysis:

# by each reaction
df.res <- get.diff.flux(imat.res1$result.model, imat.res2$result.model, nc=nc)
# total flux by each metabolite
df.met.res <- get.diff.flux.by.met(imat.res1$result.model, imat.res2$result.model, nc=nc)
# net transport flux between the cytosol and the extracellular space by each metabolite
df.tx.res <- get.diff.transport.flux(imat.res1$result.model, imat.res2$result.model, nc=nc)

df.res <- list(rxn=df.res, met=df.met.res, tx=df.tx.res)
saveRDS(df.res, file=paste0("dflux.res.",id,".RDS"))
