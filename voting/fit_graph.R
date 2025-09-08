rm(list = ls()); gc(reset = TRUE)
path <- "~/voting/"
setwd(path)
thres <- 0.0
library(abess)
library(dplyr)

load("senate_model_data.rda")
senate_nodewise <- lapply(senate_model_data, function(x) {
  dat <- t(x[["vote"]])
  splicing_gic_logistic <- slide(dat, weight = rep(1, nrow(dat)), 
                                 tune.type = "gic", ic.scale = 1, 
                                 graph.threshold = thres)
  splicing_gic_logistic[[1]]
})
save(senate_nodewise, file = "senate_nodewise.rda")
