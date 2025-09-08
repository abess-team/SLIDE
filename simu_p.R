rm(list = ls()); gc(reset = TRUE)
library(stringr)
path <- here::here(); setwd(path)
source("method_implementation.R")
source("evaluation.R")
source("simulation_main.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  type_list <- c(1)
  method <- c("SLIDE")
} else {
  for (arg in args) {
    if (grepl("^--type=", arg)) {
      type_list <- as.integer(sub("^--type=", "", arg))
    }
    if (grepl("^--method=", arg)) {
      method <- sub("^--method=", "", arg)
    }
  }
}

nrep <- 45
isparallel <- TRUE
ncore <- 5
save <- FALSE

n_max <- 1e9
if (type_list %in% c(3)) {
  # 4NN
    p_list <- c(9, 16, 25, 36, 49, 64, 81, 100)
  if (method == 'SLIDE') {
    n_start <- 1674
  } else if (method %in% c('RPLE_thres', 'RISE_thres', "logRISE_thres")) {
    n_start <- 800
  } else if (method == "ELASSO_thres") {
    n_start <- 1000
  }
  degree_list <- c(4)
} else if (type_list %in% c(1)) {
  # RRG
  p_list <- c(8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52)
  if (method %in% c('SLIDE', 'RPLE_thres', 'RISE_thres', "logRISE_thres")) {
    n_start <- 800
  } else if (method == "ELASSO_thres") {
    n_start <- 1800
  }
  degree_list <- c(3)
}
omega <- 1.5
alpha_list <- 0.3
third_order <- TRUE
beta_list <- (omega - alpha_list) / (degree_list - 1)

result_file <- "/result_p/"

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
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), 
                               t(outer(method_name, loss_name, paste0)),
                               paste0(method_name, "prop"), 
                               paste0(method_name, "mse"))
    seed_scale <- 1
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      if(isparallel == TRUE && nrep > 1) {
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, seed_scale * (1:nrep), sim, 
                          type = type, n = n_temp, p = p, method = method, 
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
        n_start_inner <- n_temp + 1
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      n_increase_scale <- 1
      if (method == "ELASSO_thres") {
        n_increase_scale <- 4
      }
      if (type == 3) {
        n_temp <- n_temp + max(round(n_increase_scale * 65536 / (n_temp + 1)), 1)
      } else if (type == 1) {
        n_temp <- n_temp + max(round(n_increase_scale * 16384 / (n_temp + 1)), 1)
      }
      
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
