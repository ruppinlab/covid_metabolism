source("./functions.R")
nc <- 4L # number of CPUs/cores

mta.res <- readRDS("../GEM/mta.res.RDS")
tmp <- readRDS("../GEM_remdesivir/mta.res.RDS")
mta.res$remdesivir <- tmp

# top 10% mta predictions
top10.mta.hits <- lapply(mta.res, get.top.hits)
# top 20% mta predictions
top20.mta.hits <- lapply(mta.res, get.top.hits, tp=20)

# enrichment of top 10% mta hits with various validation sets on the reaction level
enr.res10.rxn <- mclapply(mta.res, function(dat) lapply(vrxns, function(gs) run.enr(dat, gs, tp=10, nc=nc)), mc.cores=nc)
# enrichment of top 20% mta hits with various validation sets on the reaction level
enr.res20.rxn <- mclapply(mta.res, function(dat) lapply(vrxns, function(gs) run.enr(dat, gs, tp=20, nc=nc)), mc.cores=nc)


# AUROC, AUPRC based on the genetic screen data
auc.res <- lapply(mta.res, get.auc)

# pathway enrichment of top 10% mta hits
path.enr.res10 <- mclapply(mta.res, run.enr, gs=rset, tp=10, nc=nc, mc.cores=nc)
# pathway enrichment of top 20% mta hits
path.enr.res20 <- mclapply(mta.res, run.enr, gs=rset, tp=20, nc=nc, mc.cores=nc)

save(top10.mta.hits, top20.mta.hits, enr.res10.rxn, enr.res20.rxn, auc.res, path.enr.res10, path.enr.res20, file="check.mta.res.RData")


