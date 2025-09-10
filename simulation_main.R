# ==============================================================================
# MAIN SIMULATION FUNCTION FOR GRAPHICAL MODEL COMPARISON
# ==============================================================================
# This file implements the main simulation framework for comparing different
# graphical model estimation methods on synthetic data with known ground truth.

#' On Simulation Study for Graphical Model Methods
#' 
#' Runs a complete simulation comparing multiple graphical model estimation methods
#' including ISO variants (RPLE, RISE, logRISE), ELASSO, LogRelax, and SLIDE.
#' Generates synthetic data with known ground truth and evaluates method performance.
#' 
#' @param seed Random seed for reproducibility
#' @param type Graph structure type for data generation
#' @param n Sample size for training and validation sets
#' @param p Number of variables (nodes in the graph)
#' @param method Vector of method names to compare
#' @param alpha Signal strength parameter for data generation (default: 0.4)
#' @param beta Additional parameter for certain graph types (default: NULL)
#' @param degree Average node degree for random graphs (default: 3)
#' @param thres Global threshold for sparsity control (default: NULL, uses alpha/2)
#' @return Named vector containing all performance metrics for all methods
sim <- function(seed, type, n, p, method, alpha = 0.4, beta = NULL, degree = 3, thres = NULL) {
  runtime <- 0.0  # Track computation time for each method
  print(paste0("seed: ", seed))
  
  # DATA GENERATION
  # Generate synthetic graphical model data with known ground truth
  if (p <= 20) {
    # For small problems: use exact sampling
    train <- generate.bmn.data(n, p, type = type, seed = seed, graph.seed = 1, 
                               alpha = alpha, beta = beta, degree = degree)
    valid <- generate.bmn.data(n, theta = train[["theta"]], seed = 10000 + seed)
  } else {
    # For large problems: use Gibbs sampling for memory efficiency
    train <- generate.bmn.data(n, p, type = type, seed = seed, graph.seed = 1, 
                               alpha = alpha, beta = beta, degree = degree, method = "gibbs")
    valid <- generate.bmn.data(n, theta = train[["theta"]], seed = 10000 + seed, method = "gibbs")
  }
  
  # Extract components from generated data
  sample_weight <- c(train[["weight"]], valid[["weight"]])  # Observation frequency weights
  theta <- train[["theta"]]              # True interaction matrix (ground truth)
  edge_num_vec <- rowSums(theta != 0)    # Number of edges per node
  train <- train[["data"]]               # Training data matrix
  valid <- valid[["data"]]               # Validation data matrix
  
  # Display ground truth statistics
  true_size <- 0.5 * sum(theta != 0)  # Number of edges (accounting for symmetry)
  print(paste0("true size: ", true_size))
  
  # REGULARIZATION PARAMETER SETUP
  # Create geometric sequence of regularization parameters
  cc <- 12
  shrink <- 0.8
  c_vector <- c()
  while (cc >= 1) {
    c_vector <- c(c_vector, cc)
    cc <- floor(cc * shrink)
  }
  
  # Set default threshold if not provided
  if (is.null(thres)) {
    thres <- alpha / 2
  }
  
  # Method-specific regularization parameters (tuned based on the recommendation from Lokhov, Andrey Y., et al (2018))
  c_RPLE <- 0.2      # Regularization for RPLE
  c_RISE <- 0.4      # Regularization for RISE  
  c_logRISE <- 0.8   # Regularization for logRISE
  
  # Comprehensive grid for cross-validation
  c_list <- 10^(seq(-4, 1, length.out = min(c(p - 2, 100))))
  
  # Prepare data format: add weights as first column
  train <- cbind(sample_weight[1:nrow(train)], train)
  valid <- cbind(sample_weight[(nrow(train) + 1):(nrow(train) + nrow(valid))], valid)
  
  # ===========================================================================
  # METHOD EXECUTION SECTION
  # Each method is run conditionally based on the 'method' parameter
  # Runtime is tracked for each method execution
  # ===========================================================================
  
  # RPLE METHOD (Regularized Pseudo-Likelihood Estimator)
  if ("RPLE" %in% method)
    runtime <- system.time(
      RPLE <- ISO(train, valid, c_list = c_RPLE, method = "RPLE")
    )[3]
  else
    RPLE <- NULL
    
  # RPLE with manual thresholding
  if ("RPLE_thres" %in% method)
    runtime <- system.time(
      RPLE_thres <- ISO(train, valid, c_list = c_RPLE, thres = thres, method = "RPLE")
    )[3]
  else
    RPLE_thres <- NULL
    
  # RPLE with cross-validation and thresholding
  if ("RPLE_cv_thres" %in% method)
    runtime <- system.time(
      RPLE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "RPLE")
    )[3]
  else
    RPLE_cv_thres <- NULL
  
  # ELASSO METHOD (Extended LASSO with EBIC selection)
  if ("ELASSO" %in% method)
    runtime <- system.time(ELASSO <- ELASSO(train, valid))[3]
  else
    ELASSO <- NULL
    
  # ELASSO with manual thresholding
  if ("ELASSO_thres" %in% method)
    runtime <- system.time(
      ELASSO_thres <- ELASSO(train, valid, thres = thres)
    )[3]
  else
    ELASSO_thres <- NULL
  
  # RISE METHOD (Regularized Interaction Screening Estimator)
  if ("RISE" %in% method)
    runtime <- system.time(
      RISE <- ISO(train, valid, c_list = c_RISE, method = "RISE")
    )[3]
  else
    RISE <- NULL
    
  # RISE with manual thresholding
  if ("RISE_thres" %in% method)
    runtime <- system.time(
      RISE_thres <- ISO(train, valid, c_list = c_RISE, thres = thres, method = "RISE")
    )[3]
  else
    RISE_thres <- NULL
    
  # RISE with cross-validation and thresholding
  if ("RISE_cv_thres" %in% method)
    runtime <- system.time(
      RISE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "RISE")
    )[3]
  else
    RISE_cv_thres <- NULL
  
  # logRISE METHOD (Log-scale RISE for numerical stability)
  if ("logRISE" %in% method) {
    runtime <- system.time(
      logRISE <- ISO(train, valid, c_list = c_logRISE, method = "logRISE")
    )[3]
  } else {
    logRISE <- NULL
  }
  
  # logRISE with manual thresholding
  if ("logRISE_thres" %in% method) {
    runtime <- system.time(
      logRISE_thres <- ISO(train, valid, c_list = c_logRISE, thres = thres, method = "logRISE")
    )[3]
  } else {
    logRISE_thres <- NULL
  }
  
  # logRISE with cross-validation and thresholding
  if ("logRISE_cv_thres" %in% method) {
    runtime <- system.time(
      logRISE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "logRISE")
    )[3]
  } else {
    logRISE_cv_thres <- NULL
  }
  
  # LogRelax METHOD (Logarithmic Relaxation with empirical alpha)
  if ("LogRelax" %in% method)
    runtime <- system.time(LogRelax <- LogRelax(train, valid))[3]
  else
    LogRelax <- NULL
    
  # LogRelax with theoretical alpha adjustment
  if ("LogRelaxTAlpha" %in% method)
    runtime <- system.time(LogRelaxTAlpha <- LogRelax(train, valid, theory_alpha = TRUE))[3]
  else
    LogRelaxTAlpha <- NULL
    
  # LogRelax with manual thresholding
  if ("LogRelax_thres" %in% method)
    runtime <- system.time(LogRelax_thres <- LogRelax(train, valid, thres = thres))[3]
  else
    LogRelax_thres <- NULL
    
  # LogRelax with theoretical alpha and thresholding
  if ("LogRelaxTAlpha_thres" %in% method)
    runtime <- system.time(LogRelaxTAlpha_thres <- LogRelax(train, valid, theory_alpha = TRUE, thres = thres))[3]
  else
    LogRelaxTAlpha_thres <- NULL
  
  # SLIDE METHODS (requires different data format without weights column)
  train <- train[, -1]  # Remove weight column for SLIDE
  valid <- valid[, -1]
  pool_data <- rbind(train, valid)  # Combine all data for SLIDE
  max_size_list <- rep(max(edge_num_vec), p)
  # SLIDE with oracle knowledge (knows true support sizes)
  if ("SLIDE_oracle" %in% method) {
    runtime <- system.time(
      SLIDE_oracle <- slide(pool_data, weight = sample_weight, max.support.size = edge_num_vec, tune.type = "gic", support.size = edge_num_vec)[[1]]
    )[3]
    diag(SLIDE_oracle) <- 0.0  # Remove diagonal elements
  } else {
    SLIDE_oracle <- NULL
  }
  
  # SLIDE with automatic model selection
  if ("SLIDE" %in% method) {
    runtime <- system.time(
      SLIDE <- slide(pool_data, weight = sample_weight, tune.type = "gic", ic.scale = 1.0, graph.threshold = thres, max.support.size = max_size_list)[[1]]
    )[3]
    diag(SLIDE) <- 0.0  # Remove diagonal elements
  } else {
    SLIDE <- NULL
  }
  
  # ===========================================================================
  # RESULTS COMPILATION AND EVALUATION
  # ===========================================================================
  
  # Collect all method results in a named list
  res <- list(
    RPLE = RPLE,
    RPLE_thres = RPLE_thres,
    RPLE_cv_thres = RPLE_cv_thres,
    ELASSO = ELASSO,
    ELASSO_thres = ELASSO_thres,
    RISE = RISE,
    RISE_thres = RISE_thres,
    RISE_cv_thres = RISE_cv_thres,
    logRISE = logRISE,
    logRISE_thres = logRISE_thres,
    logRISE_cv_thres = logRISE_cv_thres,
    LogRelax = LogRelax,
    LogRelaxTAlpha = LogRelaxTAlpha,
    LogRelax_thres = LogRelax_thres,
    LogRelaxTAlpha_thres = LogRelaxTAlpha_thres,
    SLIDE_oracle = SLIDE_oracle,
    SLIDE = SLIDE
  )
  
  # PERFORMANCE EVALUATION
  # Apply all evaluation metrics to compare methods
  
  # Structure recovery performance (TPR, FPR, MCC)
  rate <- sim.rate(res, theta)
  print(rate)
  
  # Matrix norm losses (operator, L1, Frobenius)
  loss <- sim.loss(res, theta)
  
  # Exact recovery indicators (binary success/failure)
  is.rec <- sim.prop(res, theta)
  
  # Mean squared error
  mse <- sim.mse(res, theta)
  
  # Conditional AUC computation (only when no hard thresholding)
  if (thres == 0.0) {
    auc <- sim.auc(res, theta)
    res <- c(rate, loss, is.rec, mse, auc, runtime)
  } else {
    res <- c(rate, loss, is.rec, mse)
  }
  
  return(res)
}
