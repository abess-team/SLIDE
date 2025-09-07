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
  S[rank(abs(S), ties.method = "first") <= (length(S) - size)] <- 0
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
    S <- t(conf_all_centered) %*% diag(freq_vec_all) %*% conf_all_centered
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
  
  theta_est <- - solve(W)
  diag(theta_est) <- mu_bar
  if (!magnetic) {
    diag(theta_est) <- 0.0
  }
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
  return(theta_est)
}

