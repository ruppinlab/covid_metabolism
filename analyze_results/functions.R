library(my.utils)
library(gembox)
data("recon1")

load("../data/data.for.validation.RData")

# function for load mta result
load.mta.res <- function(fn, x="result.rmetal", m=recon1) {
  if (endsWith(fn, "RData")) {
    load(fn)
  } else mta.res <- readRDS(fn)
  res <- mta.res[[x]]
  res[, rank.rmta:=rank(-rTS, na.last="keep")]
  res[, norm.score:=rTS-rTS[rxn=="ctrl"]]
  res[, flag:=norm.score>0]
  res <- res[rxn!="ctrl"]
  res <- res[, .(reaction=rxn, score=rTS, norm.score, top.percent=100*rank.rmta/(.N), flag, genes=rxns2genes(m,id), pathway=m$subSystems[id])][order(top.percent)]
}

# function to get top mta hits top tp percent
get.top.hits <- function(dat, tp=10) {
  if (is.null(dat)) return(NULL)
  dat[top.percent<tp & flag==TRUE, -"flag"]
}

# pathway of reactions in recon1
rset <- subsystems2gsets(recon1)

# function to run enrichment test on mta result
run.enr <- function(dat, gs, tp=10, nc=4L, on.gene=FALSE, m=recon1) {
  if (is.null(dat)) return(NULL)
  dt <- get.top.hits(dat, tp)
  if (on.gene) {
    res <- enrich.gsets(unique(unlist(dt$genes)), gs, m$genes, nc=nc)
  } else res <- enrich.gsets(dt$reaction, gs, dat$reaction, nc=nc, name="reaction")
  res
}

# function to run gsea on top mta result against ranked gene list from genetic screen results
run.gsea1 <- function(dats.mta, tp=10, dat.rnk=screen.full) {
  dats.mta <- dats.mta[!sapply(dats.mta, is.null)]
  gs <- lapply(dats.mta, function(x) {
    dt <- get.top.hits(x, tp)
    unique(unlist(dt$genes))
  })
  mclapply(dat.rnk, function(xx) gsea(xx, gs), mc.cores=3L)
}

# function to run gsea on mta result
run.gsea <- function(dat, gs) {
  if (is.null(dat)) return(NULL)
  gsea(dat, gs, x="norm.score", id="reaction")
}

# function to get AUROC, AUPRC
get.auc <- function(dat) {
  if (is.null(dat)) return(NULL)
  tmp <- cn("crispr.wei","crispr.daniloski","rnai.sumit","union","crispr.union","rnai.sumit1","union1","rnai.sumit2","union2")
  lapply(tmp, function(x) {
    if (x %in% c("union","crispr.union","union1","union2")) {
      a <- paste0("pos.", x)
      b <- paste0("neg.", x)
    } else {
      a <- paste0(x, ".pos")
      b <- paste0(x, ".neg")
    }
    data.table(auroc=get.roc(dat[!is.na(score),score], pos=rxns.screen.for.roc[[a]], neg=rxns.screen.for.roc[[b]], x.names=dat[!is.na(score),reaction], curve=FALSE),
               auprc=get.prc(dat[!is.na(score),score], pos=rxns.screen.for.roc[[a]], neg=rxns.screen.for.roc[[b]], x.names=dat[!is.na(score),reaction], curve=FALSE))
  })
}


