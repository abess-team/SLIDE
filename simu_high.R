rm(list = ls()); gc(reset = TRUE)
args <- commandArgs(trailingOnly = TRUE)
# path <- "D:/huatong/Ising202009/code/split-code/large_scale/"
# path <- "/public/home/zhujin/Ising_xuanyu/large_scale/"
path <- "/Users/zhujin/splicing-ising/code-simulate/large-scale"
# path <- "D:/ising-L0"
setwd(path)
dir_name <- "result_high"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

source("method_implementation.R")
source("evaluation.R")
source("simulation_main.R")

save <- FALSE
# method <- c(
#   "RISE_cv_thres",
#   "logRISE_cv_thres",
#   "nodewise_logistic_gic2",
#   "ELASSO_thres",
#   "LogRelax"
# )
method <- args[1]
type_list <- c(8, 9, 10, 11, 12)
# type_list <- c(14)
## dont run type = 13 because the dimension is not match
# p_list <- c(16, 24, 32, 40, 48)
# p_list <- c(12, 16, 20, 24, 28, 32, 36, 40)
p_list <- c(18, 20, 22, 24, 26, 28, 30, 32, 34)
# p_list <- c(100, 144, 196)
alpha <- c(0.4)
beta <- c(0.5)
degree <- 3
n <- 200

method_name <- paste0(method, ".")
res_num <- length(method)
rate_name <- c("TPR", "FPR", "MCC")
loss_name <- c("loss_op", "loss_l1", "loss_F")
mark <- "test"


isparallel <- FALSE
nrep <- 48
ncore <- 6
if (isparallel) {
  library(parallel)
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  clusterExport(cl, ls())
  suppressMessages(clusterEvalQ(cl, expr = {
    library(glmnet)
    library(ROI)
    library(abess)
    library(CVXR)
  }))
} else {
  library(glmnet)
  library(ROI)
  library(abess)
}

for (p in p_list) {
  time_flag <- proc.time()[3][[1]]
  res_summary <- matrix(0, nrow = length(type_list), ncol = 10 * res_num)
  colnames(res_summary) <- c(
    t(outer(method_name, rate_name, paste0)),
    t(outer(method_name, loss_name, paste0)),
    paste0(method_name, "prop"),
    paste0(method_name, "mse"),
    paste0(method_name, "auc"),
    "runtime"
  )
  res <- list()
  for (k in 1:length(type_list)) {
    type <- type_list[k]
    print(paste0("p = ", p, "; type = ", type))
    if (isparallel) {
      clusterExport(cl,
                    c("type", "n", "p", "method", "alpha", "beta", "degree"))
      res[[k]] <- parSapply(
        cl,
        1:nrep,
        sim,
        type = type,
        n = n,
        p = p,
        method = method,
        alpha = alpha,
        beta = beta,
        degree = degree,
        thres = 0.0
      )
    } else {
      res[[k]] <- sapply(
        1,
        sim,
        type = type,
        n = n,
        p = p,
        method = method,
        alpha = alpha,
        beta = beta,
        degree = degree,
        thres = 0.0
      )
    }
    res_summary[k, ] <- rowMeans(res[[k]], na.rm = TRUE)
  }
  print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))
  write.csv(t(res_summary),
            paste0(path, result_file, method, "_p", p, "_", mark, ".csv"))
}

if (isparallel) {
  stopCluster(cl)
}

