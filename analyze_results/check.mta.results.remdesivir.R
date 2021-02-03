source("./functions.R")
nc <- 4L # number of CPUs/cores

mta.res <- readRDS("../GEM_remdesivir/mta.res.RDS")

# top 10% mta predictions
top10.mta.hits <- get.top.hits(mta.res)
# top 20% mta predictions
top20.mta.hits <- get.top.hits(mta.res, tp=20)

# enrichment of top 10% mta hits with various validation sets on the reaction level
enr.res10.rxn <- lapply(vrxns, function(gs) run.enr(mta.res, gs, tp=10, nc=nc))
# enrichment of top 20% mta hits with various validation sets on the reaction level
enr.res20.rxn <- lapply(vrxns, function(gs) run.enr(mta.res, gs, tp=20, nc=nc))

# AUROC, AUPRC based on the genetic screen data
auc.res <- get.auc(mta.res)

# pathway enrichment of top 10% mta hits
path.enr.res10 <- run.enr(mta.res, gs=rset, tp=10, nc=nc, mc.cores=nc)
# pathway enrichment of top 20% mta hits
path.enr.res20 <- run.enr(mta.res, gs=rset, tp=20, nc=nc, mc.cores=nc)

save(top10.mta.hits, top20.mta.hits, enr.res10.rxn, enr.res20.rxn, auc.res, path.enr.res10, path.enr.res20, file="check.mta.res.remdesivir.RData")


