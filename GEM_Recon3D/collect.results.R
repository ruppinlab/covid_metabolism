library(gembox)
library(stringr)
load("../data/Recon3D.RData")

# collect the differential flux analysis and MTA results for each dataset
# save the results from all datasets into a single data file for the convenience of downstream analysis

### 1. differential flux results

fns <- dir(pattern="^dflux\\.res")
names(fns) <- str_sub(fns, 11, -5)
tmp <- lapply(fns, readRDS)
df.res <- lapply(tmp, function(x) x$rxn)
df.met.res <- lapply(tmp, function(x) x$met)
df.tx.res <- lapply(tmp, function(x) x$tx)
save(df.res, df.met.res, df.tx.res, file="dflux.res.RData")

### 2. MTA results

fns <- dir(pattern="^mta\\.res")
names(fns) <- str_sub(fns, 9, -5)

load.mta.res <- function(fn, x="result.rmta", s="rTS", m=recon3d) {
  # function for load mta result
  mta.res <- readRDS(fn)
  res <- mta.res[[x]]
  res[, y:=res[[s]]]
  res[, rank.rmta:=rank(-y, na.last="keep")]
  res[, norm.score:=y-y[rxn=="ctrl"]]
  res[, flag:=norm.score>0]
  res <- res[rxn!="ctrl"]
  res <- res[, .(reaction=rxn, score=y, norm.score, top.percent=100*rank.rmta/(.N), flag, genes=rxns2genes(m,id), pathway=m$subSystems[id])][order(top.percent)]
}

mta.res <- lapply(fns, load.mta.res)
saveRDS(mta.res, file="mta.res.RDS")

