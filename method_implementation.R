library(glmnet)
library(ROI)
library(abess)
library(CVXR)

get_chisq_stat <- function(conf_mat, freq_vec) {
  p <- ncol(conf_mat)
  n <- sum(freq_vec)
  chisq_matrix <- matrix(0, p, p)
  chisq_res <- c()
  s <- 1
  
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      tab <- matrix(0, 2, 2)
      tab[1, 1] <- sum(freq_vec[conf_mat[, i] == - 1 & conf_mat[, j] == - 1])
      tab[1, 2] <- sum(freq_vec[conf_mat[, i] == - 1 & conf_mat[, j] == 1])
      tab[2, 1] <- sum(freq_vec[conf_mat[, i] == 1 & conf_mat[, j] == - 1])
      tab[2, 2] <- sum(freq_vec[conf_mat[, i] == 1 & conf_mat[, j] == 1])
      temp <- chisq.test(tab)
      chisq_res[s] <- temp[["p.value"]]
      chisq_matrix[i, j] <- temp[["statistic"]]
      s <- s + 1
    }
  }
  chisq_res <- p.adjust(chisq_res, method = "fdr")
  chisq_matrix <- chisq_matrix / n
  return(list(chisq_matrix = chisq_matrix, off_diag_num = sum(chisq_res < 0.05)))
}

get_act_set <- function(S, size) {
  # S_vector <- as.vector(abs(S))
  # cut_point <- sort(S_vector, decreasing = TRUE)[size + 1]
  # S[abs(S) <= cut_point] <- 0
  
  S[rank(abs(S), ties.method = "first") <= (length(S) - size)] <- 0
  
  # SS <- S
  # SS[SS != 0 ] <- 1
  # print(SS)
  
  act_set <- which(S != 0, arr.ind = TRUE)
  act_set <- act_set - 1
  
  return(act_set)
}

thres_theta <- function(theta, thres = NULL, cluster = NULL) {
  if(!is.null(thres)) {
    theta[abs(theta) <= thres] <- 0
  }
  p <- ncol(theta)
  
  if(!is.null(cluster)) {
    theta_vec <- as.vector(theta)
    model <- kmeans(theta_vec, cluster, iter.max = 20, nstart = 10)
    cl <- which.min(abs(model$centers))
    theta_vec[model$cluster == cl] <- 0
    theta <- matrix(theta_vec, p, p)
  }
  return(theta)
}

# opt_method: 1: line_search = false, 2: line search = true, 3: select the stepsize by hessian 
# c_type == 1: change immediately, 2: find the act_set that maximize the descent
splicing_whole <- function(train, valid, lambda = 4, gamma = 0.0, type = NULL, 
                           cluster = NULL, thres = NULL, cv = FALSE, c_vector, c_type = 1, max_iter = 1000, 
                           warmstart = TRUE, true_size = NULL, seed = NULL, 
                           opt_method = 3, inner_alpha = 0.1, inner_beta = 0.5, inner_warmstart = TRUE, save = FALSE) {
  time_flag <- proc.time()[3][[1]]
  conf_train <- train[, - 1]
  conf_valid <- valid[, - 1]
  
  freq_vec_valid <- valid[, 1]
  n_valid <- sum(freq_vec_valid)
  freq_vec_train <- train[, 1]
  n_train <- sum(freq_vec_train)
  
  p <- ncol(conf_train)
  S_res <- get_chisq_stat(conf_train, freq_vec_train)
  S <- S_res[["chisq_matrix"]]
  # print(paste0("chi-square matrix: "))
  # print(S)
  off_diag_num <- S_res[["off_diag_num"]]
  est_size <- off_diag_num
  
  if(!is.null(true_size)) {
    # s_list <- seq(from = floor(true_size - p^(1/3)), to = ceiling(true_size + p^(1/3)), by = 1)    
    s_list <- round(seq(from = true_size - 2, to = true_size + 4, by = 1))
  } else {
    len <- round(4 * sqrt(p))
    # gap = 0.5 * sqrt(p)
    if(type == 1 | type == 5) {
      s_list <- round(seq(from = 1, to = 2 * p, length.out = len))
    } else if(type == 3) {
      s_list <- round(seq(from = round(0.8 * p), to = 3 * p, length.out = len))
    } else if(type == 4) {
      s_list <- round(seq(from = 2 * p, to = round(4.3 * p), length.out = len))
    }
  }
  if(!cv) {s_list <- true_size}
  
  print(paste0("size list: ", paste0(s_list, collapse = " ")))
  
  valid_loss <- c()
  theta_hat <- list()
  theta_hat[[1]] <- matrix(0, p, p)
  
  if(!cv) {
    # eps <- 1e-2 / (n_train ^ (1 / 2))
    eps <- min(0.1 / n_train, 1e-5)
  } else {
    eps <- min(0.1 / n_train, 1e-5)
  }
  
  inner_eps <- 0.1 * eps
  # temp <- 10 * (beta - 0.5)
  # lambda <- min(4 * sqrt(2) ^ temp, 32)
  print(paste0("eps: ", eps, "   inner_eps: ", inner_eps, "   lambda: ", lambda))
  
  if(length(s_list) > 1) {
    valid_loss <- c()
    for(i in 1: length(s_list)) {
      s <- s_list[i]
      print(paste0("current size: ", s))
      act_set <- get_act_set(S, s)
      
      # 这个warmstart不太有用
      # if(warmstart) {ini_mat <- theta_hat[[i]]} else {ini_mat <- theta_hat[[1]]}
      ini_mat <- theta_hat[[1]]
      
      theta_hat[[i + 1]] <- IsingL0_inner(conf_train, freq_vec_train, n_train, ini_mat, act_set,
                                          c_vector = c_vector, c_type = c_type, opt_method = opt_method, 
                                          inner_warmstart = inner_warmstart, eps = eps, inner_eps = inner_eps, 
                                          lambda = lambda, gamma = gamma, alpha = inner_alpha, beta = inner_beta)
      
      exp_odd <- compute_exp_odd(conf_valid, theta_hat[[i + 1]])
      valid_loss[i] <- - compute_log_likelihood(conf_valid, freq_vec_valid, n_valid,  
                                                theta_hat[[i + 1]], exp_odd, gamma = 0)
      print(paste0("current valid loss: ", valid_loss[i]))
    }
    names(valid_loss) <- s_list
    print(paste0("valid loss list: ", paste0(valid_loss, collapse = " ")))
    ind <- which.min(valid_loss)
    print(valid_loss - valid_loss[ind], digits = 14)
    print("minimal difference: ")
    print(min(valid_loss[- ind] - valid_loss[ind]), digits = 14)
  } else {
    ind <- 1
  }
  
  s.opt <- s_list[ind]
  print(paste0("optimal size: ", s.opt))
  act_set <- get_act_set(S, s.opt)
  
  # act_set <- as.matrix(read.csv("out.csv", header = FALSE))
  # print("initial act_set: ")
  # print(act_set)
  
  # if(!warmstart | length(s_list) == 1) ind <- 0
  if(length(s_list) == 1) ind <- 0
  freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  conf_mat <- rbind(conf_train, conf_valid)
  n_all <- n_train + n_valid
  theta <- IsingL0_inner(conf_mat, freq_vec_all, n_all, theta_hat[[ind + 1]], act_set,
                         c_vector = c_vector, c_type = c_type, opt_method = opt_method, 
                         inner_warmstart = inner_warmstart, eps = eps, inner_eps = inner_eps, 
                         lambda = lambda, gamma = gamma, alpha = inner_alpha, beta = inner_beta)
  
  exp_odd <- compute_exp_odd(conf_mat, theta)
  final_loss <- - compute_log_likelihood(conf_mat, freq_vec_all, n_all, theta, exp_odd, gamma = 0)
  print(paste0("final total loss: ", final_loss))
  
  theta_est <- thres_theta(theta = theta, thres = thres, cluster = cluster)
  
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_train, "_", mark, "_theta_est.csv"))
  }
  print(paste0("time for splicing:", proc.time()[3][[1]] - time_flag, ' s'))
  return(theta = theta_est)
}

sym_theta <- function(theta, method) {
  p <- ncol(theta)
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(method == "max") {
        tmp <- max(theta[i, j], theta[j, i])
      } else if(method == "min") {
        tmp <- min(theta[i, j], theta[j, i])
      } else if(method == "average") {
        tmp <- (theta[i, j] + theta[j, i]) / 2
      }
      theta[i, j] <- tmp
      theta[j, i] <- tmp
    }
  }
  return(theta)
}

ISO <- function(train, valid, c_list = 0.2, thres = NULL, cluster = NULL, 
                method = c("RISE", "RPLE", "logRISE"), 
                sym = TRUE, save = FALSE, magnetic = FALSE, ebic = FALSE) {
  ## remove useless observations:
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]

  conf_train <- train[, -1]
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]
  n_valid <- sum(freq_vec_valid)
  freq_vec_train <- train[, 1]
  n_train <- sum(freq_vec_train)
  p <- ncol(train) - 1
  get_lambda <- function(c, n) c * sqrt(log((p ^ 2) / 0.05) / n)

  # Minimize objective
  if(method == "RISE") {
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop((exp(- drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
        lambda * sum(abs(thetai))
    }
  } else if(method == "logRISE") {
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop(log((exp(- drop(X %*% thetai) * Y)) %*% freq_vec / n)) + 
        lambda * sum(abs(thetai))
    }
  } else if(method == "RPLE") {
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop(log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
        lambda * sum(abs(thetai))
    }
  }
  
  wrapper <- function(thetai) obj(thetai, X, Y, freq_vec, n, lambda)
  
  if(length(c_list) > 1) {
    lambda_list <- get_lambda(c_list, n_train)
    print(paste0("lambda list: ", paste0(lambda_list, collapse = " ")))
    
    min_vloss_ind <- rep(0, p)
    for(i in 1: p) {
      print(paste0("node: ", i))
      
      X <- conf_train[, - i]
      if (magnetic) {
        X <- cbind(1, X)
      }
      var_num <- ncol(X)
      Y <- conf_train[, i]
      freq_vec <- freq_vec_train
      n <- n_train
      
      if (min(table(Y)) > 1) {
        valid_loss <- rep(0, length(lambda_list))
        for(ind_lam in 1: length(lambda_list)) {
          lambda <- lambda_list[ind_lam]
          
          # if(method == "RPLE") {
            # X_temp <- X
            # if (magnetic) {
            #   X_temp <- X_temp[, -1]
            # }
            # model_node <- glmnet(x = X_temp, y = Y, family = "binomial", weights = freq_vec, 
            #                      alpha = 1, lambda = lambda, intercept = FALSE, thresh = 1e-9)
            # thetai_est <- as.vector(model_node$beta) / 2
            # if (magnetic) {
            #   thetai_est <- as.vector(coef(model_node)) / 2
            # }
          # } else if (method == "RPLE" | method == "RISE" | method == "logRISE")
          if (method == "RPLE" | method == "RISE" | method == "logRISE") 
          {
            objective <- OP(F_objective(wrapper, n =  var_num),
                            bounds = V_bound(ld = -Inf, nobj = var_num))
            solve_temp <- tryCatch({
              ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
                        control = list(eval.max = 1000, iter.max = 1000))$solution
            }, error = function(err) {
              rep(0, var_num)
            })
            solve_temp[is.na(solve_temp)] <- 0.0
            if (thres == 0.0 && lambda > 0.0) {
              solve_tmp_abs <- abs(solve_temp)
              solve_temp[solve_tmp_abs < 0.05 * max(solve_tmp_abs)] <- 0
            }
            thetai_est <- solve_temp
          }
          if (magnetic) {
            valid_loss[ind_lam] <- obj(thetai = thetai_est, 
                                       X = cbind(1, conf_valid[, - i]), Y = conf_valid[, i], 
                                       freq_vec = freq_vec_valid, n = n_valid, lambda = 0)
          } else {
            valid_loss[ind_lam] <- obj(thetai = thetai_est, 
                                       X = conf_valid[, - i], Y = conf_valid[, i], 
                                       freq_vec = freq_vec_valid, n = n_valid, lambda = 0)
          }
        }
        min_vloss_ind[i] <- which.min(valid_loss)
        print(paste0("valid loss: ", paste0(round(valid_loss, 6), collapse = " ")))
      } else {
        min_vloss_ind[i] <- length(c_list)
      }
    }
    print(paste0("indices minimize valid loss: ", paste0(min_vloss_ind, collapse = " ")))
    c_select <- c_list[min_vloss_ind]
  } else {
    c_select <- rep(c_list, p)
  }
  print(paste0("c selected: ", paste0(c_select, collapse = " ")))
  
  conf_all <- rbind(conf_train, conf_valid)
  freq_vec <- freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  n <- n_all <- n_train + n_valid
  lambda_select <- get_lambda(c_select, n)
  
  theta_est <- matrix(0, p, p)
  for(i in 1:p) {
    lambda <- lambda_select[i]
    print(paste0("node: ", i, "  lambda: ", lambda))
    
    X <- conf_all[, - i]
    if (magnetic) {
      X <- cbind(1, X)
    }
    var_num <- ncol(X)
    Y <- conf_all[, i]
    
    if (min(table(Y)) > 1) {
      # if(method == "RPLE") {
        # X_temp <- X
        # if (magnetic) {
        #   X_temp <- X_temp[, -1]
        # }
        # model_node <- glmnet(x = X_temp, y = Y, family = "binomial", weights = freq_vec, 
        #                      alpha = 1, lambda = lambda, intercept = FALSE, thresh = 1e-9)
        # theta_est[i, - i] <- as.vector(model_node$beta) / 2
        # if (magnetic) {
        #   theta_est[i, i] <- as.vector(coef(model_node))[1] / 2
        # }
      # } else if (method == "RPLE" | method == "RISE" | method == "logRISE") {
      if (method == "RPLE" | method == "RISE" | method == "logRISE") {
        objective <- OP(F_objective(wrapper, n =  var_num),
                        bounds = V_bound(ld = -Inf, nobj = var_num))
        solve_temp <- tryCatch({
          ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
                    control = list(eval.max = 1000, iter.max = 1000))$solution
        }, error = function(err) {
          rep(0, var_num)
        })
        solve_temp[is.na(solve_temp)] <- 0.0
        if (thres == 0.0 && lambda > 0.0) {
          solve_tmp_abs <- abs(solve_temp)
          solve_temp[solve_tmp_abs < 0.05 * max(solve_tmp_abs)] <- 0
        }
        if (magnetic) {
          theta_est[i, -i] <- solve_temp[-1]
          theta_est[i, i] <- solve_temp[1]
        } else {
          theta_est[i, - i] <- solve_temp
        }
      }
    } else {
      theta_est[i, - i] <- 0
    }
  }
  print(theta_est)
  print(paste0("optimal size: ", sum(theta_est != 0) / 2))
  
  if(sym) {
    theta_est <- sym_theta(theta = theta_est, method = "max")
  }
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_origin.csv"))
  }
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_est.csv"))
  }
  print(theta_est)
  return(theta_est)
}


ELASSO <- function(train, valid, c_list = c(0.002, 0.02, 0.2, 2, 20), 
                   thres = NULL, cluster = NULL, 
                   sym = TRUE, save = FALSE, magnetic = FALSE) {
  ## remove useless observations:
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]
  
  conf_train <- train[, -1]
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]
  n_valid <- sum(freq_vec_valid)
  freq_vec_train <- train[, 1]
  n_train <- sum(freq_vec_train)
  p <- ncol(train) - 1
  
  conf_all <- rbind(conf_train, conf_valid)
  freq_vec <- freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  n <- n_all <- n_train + n_valid
  
  get_lambda <- function(c, n) c * sqrt(log((p ^ 2) / 0.05) / n)
  lambda_select <- get_lambda(c_list, n)
  
  obj <- function(thetai, X, Y, freq_vec, n, lambda) {
    drop(log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
      lambda * sum(abs(thetai))
  }
  pl_obj <- function(thetai, X, Y, freq_vec, n) {
    drop(-log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec)
  }

  wrapper <- function(thetai) obj(thetai, X, Y, freq_vec, n, lambda)

  theta_est <- matrix(0, p, p)
  for(i in 1:p) {
    X <- conf_all[, - i]
    if (magnetic) {
      X <- cbind(1, X)
    }
    var_num <- ncol(X)
    Y <- conf_all[, i]
    if (min(table(Y)) > 1) {
      ebic_value <- numeric(length(lambda_select))
      thetai_list <- list()
      for (j in 1:length(lambda_select)) {
        # solve L1 optimization
        lambda <- lambda_select[j]
        objective <- OP(F_objective(wrapper, n =  var_num),
                        bounds = V_bound(ld = -Inf, nobj = var_num))
        solve_temp <- tryCatch({
          ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
                    control = list(eval.max = 1000, iter.max = 1000))$solution
        }, error = function(err) {
          rep(0, var_num)
        })
        solve_temp[is.na(solve_temp)] <- 0.0
        if (n <= 500) {
          solve_temp[abs(solve_temp) <= 1e-4] <- 0.0
        } else {
          solve_temp[abs(solve_temp) <= 1e-8] <- 0.0
        }
        thetai_list[[j]] <- solve_temp
        # compute EBIC
        gamma <- 0.25
        penalty <- log(n_all) + 2 * gamma * log(p - 1)
        if (magnetic) {
          ebic_value[j] <- pl_obj(thetai = solve_temp, 
                                  X = cbind(1, conf_all[, - i]), 
                                  Y = conf_all[, i], 
                                  freq_vec = freq_vec_all, n = n_all)
        } else {
          ebic_value[j] <- pl_obj(thetai = solve_temp, 
                                  X = conf_all[, - i], 
                                  Y = conf_all[, i], 
                                  freq_vec = freq_vec_all, n = n_all)
        }
        ebic_value[j] <- -2 * ebic_value[j] + sum(solve_temp > 0.0) * penalty
        
      }
      print(paste0("Select lambda: ", lambda_select[which.min(ebic_value)]))
      solve_temp <- thetai_list[[which.min(ebic_value)]]
      if (magnetic) {
        theta_est[i, -i] <- solve_temp[-1]
        theta_est[i, i] <- solve_temp[1]
      } else {
        theta_est[i, - i] <- solve_temp
      }
    } else {
      theta_est[i, - i] <- 0
    }
  }
  
  print(theta_est)
  print(paste0("optimal size: ", sum(theta_est != 0) / 2))
  
  if(sym) {
    theta_est <- sym_theta(theta = theta_est, method = "max")
  }
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  print(theta_est)
  return(theta_est)
}


LogRelax <- function(train, valid, alpha=0.05, theory_alpha=FALSE, 
                     thres = NULL, cluster = NULL, 
                     save = FALSE, magnetic = FALSE) {
  ## remove useless observations:
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]
  
  conf_train <- train[, -1]
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]
  freq_vec_train <- train[, 1]

  conf_all <- rbind(conf_train, conf_valid)
  freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  p <- ncol(conf_train)
  num <- sum(freq_vec_all)
  
  mu_bar <- apply(conf_all, 2, function(x) {
    sum(x * freq_vec_all) / sum(freq_vec_all)
  })
  get_lambda <- function(mu_bar, num, p, alpha, theory_alpha) {
    sigma <- sqrt(1 - mu_bar^2)
    min_sigma_prod <- outer(sigma, sigma)
    diag(min_sigma_prod) <- max(min_sigma_prod)
    min_sigma_prod <- min(min_sigma_prod)
    if (theory_alpha) {
      denominator <- sqrt(qchisq(p = 1.0 - (alpha / (2 * p)^2), df = 1))  ## theory
    } else {
      denominator <- sqrt(qchisq(p = 1.0 - alpha, df = 1))    ## empirical usage 
    }
    numerator <- min_sigma_prod * sqrt(num)
    return(denominator / numerator)
  }
  selected_lambda <- get_lambda(mu_bar, num, p, alpha, theory_alpha)
  
  compute_S <- function(conf_all, freq_vec_all, mu_bar) {
    conf_all_centered <- sweep(conf_all, 2, mu_bar)
    freq_vec_all <- freq_vec_all / sum(freq_vec_all)
    
    # the line is equal to S implemented with for-loop
    S <- t(conf_all_centered) %*% diag(freq_vec_all) %*% conf_all_centered
    
    # S <- matrix(0.0, p, p)
    # for (i in 1:nrow(conf_all_centered)) {
    #   S <- S + freq_vec_all[i] * as.matrix(conf_all_centered[i, ]) %*% t(conf_all_centered[i, ])
    # }
    return(S)
  }
  S <- compute_S(conf_all, freq_vec_all, mu_bar)
  diag(S) <- diag(S) + 1/3
  
  # Minimize objective
  ## By CVXR
  W <- Variable(p, p, PSD = TRUE)
  obj <- log_det(W)
  matrix_error_ub <- matrix(selected_lambda, p, p)
  diag(matrix_error_ub) <- 1e-6
  constr <- list(abs(S - W) <= matrix_error_ub)
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob, solver = "SCS")
  W <- result$getValue(W)
  
  ## By ROI
  # obj <- function(W) {
  #   drop(log(det(W)))
  # }
  # objective <- OP(F_objective(obj),
  #                 bounds = V_bound(ld = -Inf, nobj = var_num))
  # solve_temp <- tryCatch({
  #   ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
  #             control = list(eval.max = 1000, iter.max = 1000))$solution
  # }, error = function(err) {
  #   rep(0, var_num)
  # })
  # solve_temp[is.na(solve_temp)] <- 0.0
  
  theta_est <- - solve(W)
  diag(theta_est) <- mu_bar
  if (!magnetic) {
    diag(theta_est) <- 0.0
  }
  # theta_est[abs(theta_est) <= 1e-4] <- 0.0
  # print(round(theta_est, digits = 2))
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_origin.csv"))
  }
  
  if (is.null(thres) && selected_lambda > 0.0) {
    theta_est_abs <- abs(theta_est)
    theta_est[theta_est_abs < 0.05 * max(theta_est_abs)] <- 0.0
  } 
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_est.csv"))
  }
  # print(round(theta_est, digits = 2))
  return(theta_est)
}

