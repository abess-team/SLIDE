rm(list = ls()); gc(reset = TRUE)
path <- "/Users/zhujin/splicing-ising/code-simulate/code-github"
setwd(path)

library(stringr)
source("method_implementation.R")
source("evaluation.R")
source("simulation_main.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # method <- c("RPLE_thres")
  # method <- c("RISE_thres")
  # method <- c("logRISE_thres")
  # method <- c("ELASSO_thres")
  method <- c("SLIDE")
  type_list <- c(1)
} else {
  method <- sub("^--method=", "", args[grep("^--method=", args)])
  type_list <- as.numeric(sub("^--type=", "", args[grep("^--type=", args)]))
}

nrep <- 45
isparallel <- TRUE
ncore <- 45
save <- FALSE

p_list <- c(16)
alpha_list <- c(0.4)
degree_list <- c(3)

n_start <- 1e2
if (str_detect(method, "LogRelax")) {
  n_max <- 1e6
} else {
  n_max <- 2e9
}

if (type_list == 1) {
  # Type B
  beta_list <- (14:24)
} else if (type_list == 2) {
  # Type D
  beta_list <- (20:32)
} else if (type_list == 3) {
  # Type A
  beta_list <- (20:36) / 2
} else if (type_list == 4) {
  # Type C
  beta_list <- (22:32)
} else if (type_list == 5) {
  # Type E
  beta_list <- (20:36) / 2
}
beta_list <- beta_list / 20

dir_name <- paste0("result_type", type_list)
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

p_list <- sort(p_list)
alpha_list <- sort(alpha_list, decreasing = TRUE)
beta_list <- sort(beta_list)
degree_list <- sort(degree_list)
l1 <- length(p_list); l2 <- length(alpha_list); l3 <- length(beta_list); l4 <- length(degree_list); 
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
    library(CVXR)
  }))
}

time_flag <- proc.time()[3][[1]]
stop_flag <- FALSE
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
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      if(isparallel == TRUE && nrep > 1) {
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, 1:nrep, sim, type = type, n = n_temp, p = p, method = method, 
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
      } else {
        res_summary <- rbind(res_summary, temp)
      }
      rownames(res_summary)[k] <- paste0("n_", n_temp)
      
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      if(prop == 1) {
        n_start_inner <- n_temp
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      n_temp <- round(n_temp * (2 - prop))
      k <- k + 1
    }
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary,rep(0, 8))
      stop_flag <- TRUE
    }
    print("______________________________________________________________")
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
    if (stop_flag) {
      print("Stop sufficient sample size searching!")
      break;
    }
  }
}

if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))
