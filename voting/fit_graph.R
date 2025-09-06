rm(list = ls()); gc(reset = TRUE)
path <- "~/voting/"
setwd(path)
thres <- 0.0
library(abess)
library(dplyr)

load("senate_model_data.rda")
senate_nodewise <- lapply(senate_model_data, function(x) {
  dat <- t(x[["vote"]])
  splicing_gic_logistic <- nodewise_L0(dat, weight = rep(1, nrow(dat)), 
                                       tune.type = "gic", ic.scale = 1, 
                                       graph.threshold = thres, magnetic = TRUE)
  splicing_gic_logistic[[1]]
})
save(senate_nodewise, file = "senate_nodewise.rda")

load("senate113.rda")
dat[dat == 2] <- 0
dat[dat == 0] <- -1
splicing_gic_logistic <- nodewise_L0(dat, weight = rep(1, nrow(dat)), 
                                     tune.type = "gic", ic.scale = 1, 
                                     graph.threshold = thres, magnetic = TRUE)
splicing_gic_logistic <- splicing_gic_logistic[[1]]
save(splicing_gic_logistic, file = "res_senate113_nodewise.rda")
