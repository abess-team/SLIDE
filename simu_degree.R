## 尝试 3 比较多次项
# --> 显著，但是系数是负的
# --> 1. 增大重复次数 (check)
# --> 2. 换成不同的随机数种子 (hopefully) 
# --> 2.1 increase p to 20 (almost) but our proposal also close to quadtic...
# --> 2.2 increase gamma such that the gap of empirical sample complexity enlarges when degree increases
# --> [x] 2.2 decrease the added samples when enlarging samples (just make this more smooth...)

rm(list = ls()); gc(reset = TRUE)
args <- commandArgs(trailingOnly = TRUE)
# path <- "D:/huatong/Ising202009/code/split-code/large_scale/"
# path <- "/public/home/zhujin/Ising_xuanyu/large_scale/"
path <- "/Users/zhujin/splicing-ising/code-simulate/large-scale"
# path <- "D:/ising-L0"
setwd(path)
source("method_implementation.R")
source("evaluation.R")
source("simulation_main.R")

nrep <- 45
isparallel <- TRUE
ncore <- 5
save <- FALSE
# 只能一个 method!!!!!
# method <- c("RPLE_thres")
method <- c("ELASSO_thres")
# method <- c("nodewise_logistic_gic2")
# method <- c("nodewise_logistic_cv")
# method <- c("logRISE_thres")
# method <- args[1]

# code 里的n都是训练集大小 画图时要 * 2
n_max <- 1e9

# 除了type以外，一次只能跑一个维数 不然n_start出问题

third_order <- FALSE

##### local dense network #####
type_list <- c(18)
p_list <- c(16)
omega <- 2.0
degree_list <- c(1:14)
n_start <- 2e1
third_order <- TRUE

##### much fast (3-order or 4-order. Useless! Just ignore!) #####
# type_list <- c(16)
# p_list <- c(24)
# omega <- 1.5
# degree_list <- c(1:22)
# n_start <- 2e1
# third_order <- TRUE

##### much fast #####
# type_list <- c(16)
# p_list <- c(16)
# omega <- 2.8
# degree_list <- c(3:12)
# n_start <- 5e4
# alpha_list <- c(0.2)

##### slower #####
# type_list <- c(16)
# p_list <- c(18)
# omega <- 3.2
# degree_list <- c(3:16)
# n_start <- 5e4
# alpha_list <- c(0.2)

##### very slow #####
# type_list <- c(16)
# p_list <- c(20)
# omega <- 3.6
# degree_list <- c(4:18)
# n_start <- 1e7
# alpha_list <- c(0.2)

##### much fast #####
# type_list <- c(17)  ## spin glass is easy to hold (it may because not all of node are homogeneous)
# p_list <- c(16)
# omega <- 2.8
# degree_list <- c(3:14)
# n_start <- 5e4
# alpha_list <- c(0.2)

if (third_order) {
  alpha_list <- omega / degree_list
  beta_list <- omega / degree_list
} else {
  beta_list <- (omega - alpha_list) / (degree_list - 1)
}
result_file <- "/result_d/"

l1 <- length(p_list); l2 <- length(alpha_list); 
l3 <- length(beta_list); l4 <- length(degree_list); 
conf_mat <- cbind(p_list, alpha_list, beta_list, degree_list)

method_name <- paste0(method, ".")
res_num <- length(method)
rate_name <- c("TPR", "FPR", "MCC")
loss_name <- c("loss_op", "loss_l1", "loss_F")
mark <- "test"

if (isparallel && nrep > 1) {
  library(parallel)
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  clusterExport(cl, ls())
  suppressMessages(clusterEvalQ(cl, expr = {
    library(glmnet)
    library(ROI)
    library(abess)
  }))
}

time_flag <- proc.time()[3][[1]]

for(type in type_list) {
  print(paste0("type: ", type))
  n_start_inner <- n_start
  
  for(conf_index in 1:nrow(conf_mat)) {
    conf <- conf_mat[conf_index, ]
    p <- conf[1]; alpha <- conf[2]; beta <- conf[3]; degree <- conf[4]
    
    info <- paste0("p", p, "_alpha", alpha, "_beta", beta, "_degree", degree)
    print(paste0("config: ", info))
    
    n_temp <- n_start_inner
    prop <- 0; k <- 1
    res_summary <- matrix(0, nrow = 1, ncol = 8 * res_num)
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), t(outer(method_name, loss_name, paste0)),
                               paste0(method_name, "prop"), paste0(method_name, "mse"))
    seed_scale <- 1
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      if(isparallel == TRUE && nrep > 1) {
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, seed_scale * (1:nrep), sim, type = type, n = n_temp, p = p, method = method, 
                          alpha = alpha, beta = beta, degree = degree)
      } else {
        temp <- sapply(1, sim, type = type, n = n_temp, p = p, method = method, 
                       alpha = alpha, beta = beta, degree = degree)
      }
      print(paste0("is recovery: ", paste0(temp[7, ], collapse = " ")))
      temp <- rowMeans(temp, na.rm = TRUE)
      prop <- temp[7]
      
      if(k == 1) {
        res_summary[k, ] <- temp
        seed_scale <- 1
      } else {
        res_summary <- rbind(res_summary, temp)
        seed_scale <- seed_scale + 1
      }
      rownames(res_summary)[k] <- paste0("n_", n_temp)
      
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      if(prop == 1) {
        n_start_inner <- round(n_temp * 1.005)
        # n_start_inner <- round(n_temp * 0.4)
        # n_start_inner <- n_temp
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      n_temp <- round(n_temp * (2.0 - prop))
      k <- k + 1
    }
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary, rep(0, 8))
    }
    print("______________________________________________________________")
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
  }
}

if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))
