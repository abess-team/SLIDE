library(ROI)    # R Optimization Infrastructure
library(CVXR)   # Convex optimization
library(abess)

#' Apply Thresholding to Interaction Matrix
#' 
#' Applies hard thresholding or clustering-based thresholding to enforce
#' sparsity in the estimated interaction matrix.
#' 
#' @param theta Estimated interaction matrix
#' @param thres Hard threshold value (elements with |value| <= thres set to 0)
#' @param cluster Number of clusters for k-means thresholding (optional)
#' @return Thresholded interaction matrix
thres_theta <- function(theta, thres = NULL, cluster = NULL) {
  # Apply hard thresholding if threshold value provided
  if(!is.null(thres)) {
    theta[abs(theta) <= thres] <- 0
  }
  
  # Apply clustering-based thresholding if number of clusters specified (generally not used)
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

#' Symmetrize Interaction Matrix
#' 
#' Enforces symmetry in the interaction matrix using different aggregation methods
#' for off-diagonal elements that should be equal.
#' 
#' @param theta Asymmetric interaction matrix
#' @param method Symmetrization method: "max", "min", or "average"
#' @return Symmetric interaction matrix
sym_theta <- function(theta, method) {
  p <- ncol(theta)

  # Symmetrize off-diagonal elements
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(method == "max") {
        tmp <- max(theta[i, j], theta[j, i])    # Take maximum
      } else if(method == "min") {
        tmp <- min(theta[i, j], theta[j, i])    # Take minimum
      } else if(method == "average") {
        tmp <- (theta[i, j] + theta[j, i]) / 2  # Take average
      }
      # Assign symmetric values
      theta[i, j] <- tmp
      theta[j, i] <- tmp
    }
  }

  return(theta)
}


#' Interaction Screening Objective (ISO) for Graphical Model Estimation
#' 
#' Implements the ISO method with different loss functions (RISE, RPLE, logRISE)
#' for estimating sparse interaction matrices in Ising models.
#' 
#' @param train Training dataset (first column: weights, remaining: binary variables)
#' @param valid Validation dataset (same format as train)
#' @param c_list Vector of regularization parameter values to try
#' @param thres Hard threshold for sparsity (default: NULL)
#' @param cluster Number of clusters for k-means thresholding (default: NULL)
#' @param method Estimation method: "RISE", "RPLE", or "logRISE"
#' @param sym Whether to symmetrize the result (default: TRUE)
#' @param save Whether to save intermediate results (default: FALSE)
#' @param magnetic Whether to estimate diagonal parameters (default: FALSE)
#' @param ebic Whether to use EBIC for model selection (default: FALSE)
#' @return Estimated sparse interaction matrix
ISO <- function(train, valid, c_list = 0.2, thres = NULL, cluster = NULL, 
                method = c("RISE", "RPLE", "logRISE"), 
                sym = TRUE, save = FALSE, magnetic = FALSE, ebic = FALSE) {
  # DATA PREPROCESSING
  # Remove observations with zero weight (unobserved configurations)
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]

  # Separate configuration matrices and frequency vectors
  conf_train <- train[, -1]       # Binary samples
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]    # Observation weights/frequencies
  n_valid <- sum(freq_vec_valid)  # Total validation sample size
  freq_vec_train <- train[, 1]
  n_train <- sum(freq_vec_train)  # Total training sample size
  p <- ncol(train) - 1             # Number of variables

  # REGULARIZATION PARAMETER CALCULATION
  # Adaptive regularization based on problem dimensions and sample suggested by Lokhov, Andrey Y., et al (2018).
  get_lambda <- function(c, n) c * sqrt(log((p ^ 2) / 0.05) / n)

  # OBJECTIVE FUNCTION DEFINITIONS
  # Define different loss functions for the ISO method
  if(method == "RISE") {
    # Regularized Interaction Screening Estimator
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop((exp(- drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
        lambda * sum(abs(thetai))
    }
  } else if(method == "logRISE") {
    # Log-RISE estimator
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop(log((exp(- drop(X %*% thetai) * Y)) %*% freq_vec / n)) + 
        lambda * sum(abs(thetai))
    }
  } else if(method == "RPLE") {
    # Regularized Pseudo-Likelihood Estimator
    obj <- function(thetai, X, Y, freq_vec, n, lambda) {
      drop(log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
        lambda * sum(abs(thetai))
    }
  }
  
  # Create wrapper function for optimization
  wrapper <- function(thetai) obj(thetai, X, Y, freq_vec, n, lambda)
  
  # REGULARIZATION PARAMETER SELECTION
  if(length(c_list) > 1) {
    # Cross-validation to select optimal regularization parameter
    lambda_list <- get_lambda(c_list, n_train)
    print(paste0("lambda list: ", paste0(lambda_list, collapse = " ")))
    
    min_vloss_ind <- rep(0, p)  # Store optimal lambda index for each node
    
    # Node-wise parameter selection via validation
    for(i in 1: p) {
      print(paste0("node: ", i))
      
      # Prepare design matrix (all variables except node i)
      X <- conf_train[, -i]
      if (magnetic) {
        X <- cbind(1, X)  # Add intercept for magnetic field
      }
      var_num <- ncol(X)
      Y <- conf_train[, i]    # Response variable (node i)
      freq_vec <- freq_vec_train
      n <- n_train
      
      # Only proceed if there's variation in the response
      if (min(table(Y)) > 1) {
        valid_loss <- rep(0, length(lambda_list))
        
        # Try each regularization parameter
        for(ind_lam in 1: length(lambda_list)) {
          lambda <- lambda_list[ind_lam]
          
          # OPTIMIZATION STEP
          # Set up and solve the optimization problem
          if (method == "RPLE" | method == "RISE" | method == "logRISE") {
            objective <- OP(F_objective(wrapper, n = var_num),
                            bounds = V_bound(ld = -Inf, nobj = var_num))
            solve_temp <- tryCatch({
              ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
                        control = list(eval.max = 1000, iter.max = 1000))$solution
            }, error = function(err) {
              rep(0, var_num)  # Return zeros if optimization fails
            })
            solve_temp[is.na(solve_temp)] <- 0.0  # Handle NaN values
            
            # Apply thresholding on tiny values that comes from interaction error
            if (thres == 0.0 && lambda > 0.0) {
              solve_tmp_abs <- abs(solve_temp)
              solve_temp[solve_tmp_abs < 0.05 * max(solve_tmp_abs)] <- 0
            }
            thetai_est <- solve_temp
          }
          
          # VALIDATION STEP
          # Compute validation loss (without regularization penalty)
          if (magnetic) {
            valid_loss[ind_lam] <- obj(thetai = thetai_est, 
                                       X = cbind(1, conf_valid[, -i]), 
                                       Y = conf_valid[, i], 
                                       freq_vec = freq_vec_valid, 
                                       n = n_valid, lambda = 0)
          } else {
            valid_loss[ind_lam] <- obj(thetai = thetai_est, 
                                       X = conf_valid[, -i], 
                                       Y = conf_valid[, i], 
                                       freq_vec = freq_vec_valid, 
                                       n = n_valid, lambda = 0)
          }
        }
        
        # Select lambda with minimum validation loss
        min_vloss_ind[i] <- which.min(valid_loss)
        print(paste0("valid loss: ", paste0(round(valid_loss, 6), collapse = " ")))
      } else {
        # If no variation in response, use largest regularization
        min_vloss_ind[i] <- length(c_list)
      }
    }
    print(paste0("indices minimize valid loss: ", paste0(min_vloss_ind, collapse = " ")))
    c_select <- c_list[min_vloss_ind]  # Selected regularization parameters
  } else {
    # Use same regularization parameter for all nodes
    c_select <- rep(c_list, p)
  }
  print(paste0("c selected: ", paste0(c_select, collapse = " ")))
  
  # FINAL ESTIMATION ON COMBINED DATA
  # Combine training and validation data for final estimation
  conf_all <- rbind(conf_train, conf_valid)
  freq_vec <- freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  n <- n_all <- n_train + n_valid
  lambda_select <- get_lambda(c_select, n)  # Final regularization parameters
  
  theta_est <- matrix(0, p, p)  # Initialize interaction matrix estimate
  
  # NODE-WISE ESTIMATION
  # Estimate each row of the interaction matrix separately
  for(i in 1:p) {
    lambda <- lambda_select[i]
    print(paste0("node: ", i, "  lambda: ", lambda))
    
    # Prepare design matrix and response for node i
    X <- conf_all[, -i]
    if (magnetic) {
      X <- cbind(1, X)  # Add intercept term
    }
    var_num <- ncol(X)
    Y <- conf_all[, i]
    
    # Only estimate if there's variation in the response
    if (min(table(Y)) > 1) {
      # Solve optimization problem for node i
      if (method == "RPLE" | method == "RISE" | method == "logRISE") {
        objective <- OP(F_objective(wrapper, n = var_num),
                        bounds = V_bound(ld = -Inf, nobj = var_num))
        solve_temp <- tryCatch({
          ROI_solve(objective, solver = "nlminb", start = rep(0, var_num),
                    control = list(eval.max = 1000, iter.max = 1000))$solution
        }, error = function(err) {
          rep(0, var_num)
        })
        solve_temp[is.na(solve_temp)] <- 0.0
        
        # Apply soft thresholding if needed
        if (thres == 0.0 && lambda > 0.0) {
          solve_tmp_abs <- abs(solve_temp)
          solve_temp[solve_tmp_abs < 0.05 * max(solve_tmp_abs)] <- 0
        }
        
        # Store results in interaction matrix
        if (magnetic) {
          theta_est[i, -i] <- solve_temp[-1]  # Off-diagonal elements
          theta_est[i, i] <- solve_temp[1]    # Diagonal element (magnetic field)
        } else {
          theta_est[i, -i] <- solve_temp     # Off-diagonal elements only
        }
      }
    } else {
      # If no variation, set row to zero
      theta_est[i, -i] <- 0
    }
  }
  
  print(theta_est)
  print(paste0("optimal size: ", sum(theta_est != 0) / 2))
  
  # POST-PROCESSING
  # Apply symmetrization if requested
  if(sym) {
    theta_est <- sym_theta(theta = theta_est, method = "max")
  }
  
  # Save intermediate results if requested
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_origin.csv"))
  }
  
  # Apply thresholding for final sparsity
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  
  # Save final results if requested
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_est.csv"))
  }
  
  print(theta_est)
  return(theta_est)
}


#' Extended LASSO (ELASSO) for Graphical Model Estimation
#' 
#' Implements the Extended LASSO method using EBIC (Extended Bayesian Information Criterion)
#' for model selection in Ising graphical models.
#' 
#' @param train Training dataset (first column: weights, remaining: binary variables)
#' @param valid Validation dataset (same format as train)
#' @param c_list Vector of regularization parameter multipliers
#' @param thres Hard threshold for sparsity (default: NULL)
#' @param cluster Number of clusters for k-means thresholding (default: NULL)
#' @param sym Whether to symmetrize the result (default: TRUE)
#' @param save Whether to save intermediate results (default: FALSE)
#' @param magnetic Whether to include node-specific parameters (default: FALSE)
#' @return Estimated sparse interaction matrix using EBIC model selection
ELASSO <- function(train, valid, c_list = c(0.002, 0.02, 0.2, 2, 20), 
                   thres = NULL, cluster = NULL, 
                   sym = TRUE, save = FALSE, magnetic = FALSE) {
  # DATA PREPROCESSING
  # Remove observations with zero frequency weight
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]
  
  # Extract configuration matrices and frequency vectors
  conf_train <- train[, -1]
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]
  n_valid <- sum(freq_vec_valid)
  freq_vec_train <- train[, 1]
  n_train <- sum(freq_vec_train)
  p <- ncol(train) - 1
  
  # Combine training and validation data
  conf_all <- rbind(conf_train, conf_valid)
  freq_vec <- freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  n <- n_all <- n_train + n_valid
  
  # REGULARIZATION PARAMETER SETUP
  get_lambda <- function(c, n) c * sqrt(log((p ^ 2) / 0.05) / n)
  lambda_select <- get_lambda(c_list, n)
  
  # OBJECTIVE FUNCTION DEFINITIONS
  # Regularized objective function (with L1 penalty)
  obj <- function(thetai, X, Y, freq_vec, n, lambda) {
    drop(log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec) / n + 
      lambda * sum(abs(thetai))
  }
  # Pseudo-likelihood objective (without penalty, for EBIC computation)
  pl_obj <- function(thetai, X, Y, freq_vec, n) {
    drop(-log(1 + exp(- 2 * drop(X %*% thetai) * Y)) %*% freq_vec)
  }

  wrapper <- function(thetai) obj(thetai, X, Y, freq_vec, n, lambda)

  theta_est <- matrix(0, p, p)  # Initialize interaction matrix

  # NODE-WISE ESTIMATION WITH EBIC SELECTION
  for(i in 1:p) {
    # Prepare design matrix and response for node i
    X <- conf_all[, -i]
    if (magnetic) {
      X <- cbind(1, X)  # Add intercept for magnetic field
    }
    var_num <- ncol(X)
    Y <- conf_all[, i]

    # Only proceed if there's variation in the response
    if (min(table(Y)) > 1) {
      ebic_value <- numeric(length(lambda_select))  # Store EBIC values
      thetai_list <- list()                         # Store coefficient estimates

      # Try each regularization parameter
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

        # Apply thresholding based on sample size
        if (n <= 500) {
          solve_temp[abs(solve_temp) <= 1e-4] <- 0.0
        } else {
          solve_temp[abs(solve_temp) <= 1e-8] <- 0.0
        }
        thetai_list[[j]] <- solve_temp

        # EBIC COMPUTATION
        # Extended Bayesian Information Criterion for model selection
        gamma <- 0.25                                    # EBIC parameter
        penalty <- log(n_all) + 2 * gamma * log(p - 1)  # EBIC penalty term
        
        # Compute log-likelihood
        if (magnetic) {
          ebic_value[j] <- pl_obj(thetai = solve_temp, 
                                  X = cbind(1, conf_all[, -i]), 
                                  Y = conf_all[, i], 
                                  freq_vec = freq_vec_all, n = n_all)
        } else {
          ebic_value[j] <- pl_obj(thetai = solve_temp, 
                                  X = conf_all[, -i], 
                                  Y = conf_all[, i], 
                                  freq_vec = freq_vec_all, n = n_all)
        }

        # Compute EBIC: -2*log-likelihood + penalty*model_size
        ebic_value[j] <- -2 * ebic_value[j] + sum(solve_temp > 0.0) * penalty
        
      }

      # Select model with minimum EBIC
      print(paste0("Select lambda: ", lambda_select[which.min(ebic_value)]))
      solve_temp <- thetai_list[[which.min(ebic_value)]]

      # Store results in precision matrix
      if (magnetic) {
        theta_est[i, -i] <- solve_temp[-1]  # Off-diagonal elements
        theta_est[i, i] <- solve_temp[1]    # Diagonal element
      } else {
        theta_est[i, -i] <- solve_temp     # Off-diagonal elements only
      }
    } else {

      # If no variation in response, set row to zero
      theta_est[i, -i] <- 0
    }
  }
  
  print(theta_est)
  print(paste0("optimal size: ", sum(theta_est != 0) / 2))

  # POST-PROCESSING
  # Apply symmetrization if requested
  if(sym) {
    theta_est <- sym_theta(theta = theta_est, method = "max")
  }

  # Apply thresholding for final sparsity
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  
  print(theta_est)
  return(theta_est)
}

#' Logarithmic Relaxation (LogRelax) Method for Graphical Model Estimation
#' 
#' Implements the LogRelax method using convex relaxation of the log-determinant
#' for interaction matrix estimation in Gaussian graphical models.
#' 
#' @param train Training dataset (first column: weights, remaining: variables)
#' @param valid Validation dataset (same format as train)
#' @param alpha Significance level for constraint relaxation (default: 0.05)
#' @param theory_alpha Whether to use theoretical alpha adjustment (default: FALSE)
#' @param thres Hard threshold for sparsity (default: NULL)
#' @param cluster Number of clusters for k-means thresholding (default: NULL)
#' @param save Whether to save intermediate results (default: FALSE)
#' @param magnetic Whether to include diagonal elements (default: FALSE)
#' @return Estimated sparse interaction matrix via convex optimization
LogRelax <- function(train, valid, alpha=0.05, theory_alpha=FALSE, 
                     thres = NULL, cluster = NULL, 
                     save = FALSE, magnetic = FALSE) {
  # DATA PREPROCESSING
  # Remove observations with zero frequency weight
  train <- train[train[, 1] > 0, ]
  valid <- valid[valid[, 1] > 0, ]
  
  # Extract configuration matrices and frequency vectors
  conf_train <- train[, -1]
  conf_valid <- valid[, -1]
  freq_vec_valid <- valid[, 1]
  freq_vec_train <- train[, 1]

  # Combine training and validation data
  conf_all <- rbind(conf_train, conf_valid)
  freq_vec_all <- c(freq_vec_train, freq_vec_valid)
  p <- ncol(conf_train)        # Number of variables
  num <- sum(freq_vec_all)     # Total sample size
  
  # SAMPLE STATISTICS COMPUTATION
  # Compute weighted sample means for each variable
  mu_bar <- apply(conf_all, 2, function(x) {
    sum(x * freq_vec_all) / sum(freq_vec_all)
  })
  
  # CONSTRAINT PARAMETER CALCULATION
  # Adaptive constraint parameter based on statistical theory
  get_lambda <- function(mu_bar, num, p, alpha, theory_alpha) {
    # Compute standard deviations assuming binary variables
    sigma <- sqrt(1 - mu_bar^2)
    
    # Find minimum product of standard deviations (tightest constraint)
    min_sigma_prod <- outer(sigma, sigma)
    diag(min_sigma_prod) <- max(min_sigma_prod)  # Exclude diagonal
    min_sigma_prod <- min(min_sigma_prod)
    
    # Choose critical value based on theory or empirical usage (see Banerjee, O., et. al 2008)
    if (theory_alpha) {
      # Theoretical adjustment for multiple testing
      denominator <- sqrt(qchisq(p = 1.0 - (alpha / (2 * p)^2), df = 1))
    } else {
      # Empirical usage (more liberal)
      denominator <- sqrt(qchisq(p = 1.0 - alpha, df = 1))
    }
    
    # Scale by sample size and minimum variance product
    numerator <- min_sigma_prod * sqrt(num)
    return(denominator / numerator)
  }
  
  selected_lambda <- get_lambda(mu_bar, num, p, alpha, theory_alpha)
  
  # SAMPLE COVARIANCE COMPUTATION
  # Compute weighted sample covariance matrix
  compute_S <- function(conf_all, freq_vec_all, mu_bar) {
    # Center the data
    conf_all_centered <- sweep(conf_all, 2, mu_bar)
    # Normalize frequency weights
    freq_vec_all <- freq_vec_all / sum(freq_vec_all)
    # Compute weighted covariance
    S <- t(conf_all_centered) %*% diag(freq_vec_all) %*% conf_all_centered
    return(S)
  }
  
  S <- compute_S(conf_all, freq_vec_all, mu_bar)
  # Add small regularization to diagonal for numerical stability
  diag(S) <- diag(S) + 1/3
  
  # CONVEX OPTIMIZATION PROBLEM
  # Solve the log-determinant maximization with elementwise constraints
  W <- Variable(p, p, PSD = TRUE)  # Precision matrix variable (positive semidefinite)
  obj <- log_det(W)                # Objective: log-determinant (concave)
  
  # Constraint matrix: allow larger errors on diagonal
  matrix_error_ub <- matrix(selected_lambda, p, p)
  diag(matrix_error_ub) <- 1e-6
  
  # Elementwise constraints: |S - W| <= error bounds
  constr <- list(abs(S - W) <= matrix_error_ub)
  
  # Solve the convex optimization problem
  prob <- Problem(Maximize(obj), constr)
  result <- solve(prob, solver = "SCS")
  W <- result$getValue(W)
  
  # PRECISION MATRIX RECOVERY
  # The precision matrix is the negative inverse of W
  theta_est <- - solve(W)
  
  # Set diagonal elements based on magnetic field option
  if (magnetic) {
    diag(theta_est) <- mu_bar  # Use sample means as magnetic fields
  } else {
    diag(theta_est) <- 0.0     # No magnetic fields
  }
  
  # Save intermediate results if requested
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_origin.csv"))
  }
  
  # POST-PROCESSING
  # Apply automatic thresholding if no manual threshold specified
  if (is.null(thres) && selected_lambda > 0.0) {
    theta_est_abs <- abs(theta_est)
    # Set small elements (< 5% of maximum) to zero
    theta_est[theta_est_abs < 0.05 * max(theta_est_abs)] <- 0.0
  } 
  
  # Apply manual or clustering-based thresholding
  theta_est <- thres_theta(theta = theta_est, thres = thres, cluster = cluster)
  
  # Save final results if requested
  if(save) {
    write.csv(theta_est, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", 
                                alpha, "_beta", beta, "_n", n_all / 2, "_", mark, "_theta_est.csv"))
  }
  
  return(theta_est)
}
