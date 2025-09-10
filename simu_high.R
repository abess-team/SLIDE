# ==============================================================================
# SIMULATION STUDY: EVALUATION ON HIGH-DIMENSIONAL CASES
# ==============================================================================
# Purpose: This script evaluates graphical model estimation methods in 
# high-dimensional settings with fixed sample size but varying dimensions of  
# variables. Unlike other simulations that find minimum sample sizes,
# this tests method performance across increasing dimensions.
#
# Input: Command line arguments specifying the estimation method
# Output: CSV files containing performance metrics across dimensions and graph types
# ==============================================================================

# ENVIRONMENT INITIALIZATION
# Clear workspace and prepare clean environment
rm(list = ls()); gc(reset = TRUE)

# Set working directory and source simulation framework
path <- here::here(); setwd(path)

# Create results directory for high-dimensional experiments
dir_name <- "result_high"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

# Source all required simulation components
source("method_implementation.R")  # Implementation of estimation methods
source("evaluation.R")             # Performance evaluation functions
source("simulation_main.R")        # Main simulation control framework

# COMMAND LINE ARGUMENT PROCESSING
# Parse command line arguments to determine which method to evaluate
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Default method when no arguments provided
  # Available methods include:
  # 1. "RPLE_cv_thres"     - RPLE with cross-validation and thresholding
  # 2. "RISE_cv_thres"     - RISE with cross-validation and thresholding  
  # 3. "logRISE_cv_thres"  - logRISE with cross-validation and thresholding
  # 4. "SLIDE"             - SLIDE method
  # 5. "ELASSO_thres"      - Extended LASSO with thresholding
  # 6. "LogRelax"          - Logarithmic relaxation method
  method <- c("SLIDE")
} else {
  # Extract method name from --method= argument
  if (grepl("^--method=", arg)) {
    method <- sub("^--method=", "", arg)
  }
}

# EXPERIMENTAL CONFIGURATION
# Core simulation parameters
save <- FALSE                    # Whether to save intermediate results
type_list <- c(1, 2, 3, 4, 5)    # Multiple graph types to test comprehensively

# HIGH-DIMENSIONAL PARAMETER SWEEP
# Test increasing problem dimensions to assess scalability
p_list <- c(18, 20, 22, 24, 26, 28, 30, 32, 34)  # Number of variables (nodes)

# Fixed parameters across all high-dimensional experiments
alpha <- c(0.4)      # Minimal signal strength
beta <- c(0.5)       # Interaction strength  
degree <- 3          # Average node degree
n <- 200             # Fixed sample size (key difference from other simulations)

# RESULTS STRUCTURE CONFIGURATION
# Set up naming conventions for results organization
method_name <- paste0(method, ".")           # Method name with separator
res_num <- length(method)                    # Number of methods (should be 1)
rate_name <- c("TPR", "FPR", "MCC")         # Binary classification metrics
loss_name <- c("loss_op", "loss_l1", "loss_F")  # Matrix reconstruction losses
mark <- "test"                               # Experiment identifier

# PARALLEL PROCESSING SETUP
# Configure parallel execution for computational efficiency
isparallel <- TRUE   # Enable parallel processing
nrep <- 50          # Number of Monte Carlo replications
ncore <- 50         # Number of CPU cores to utilize

if (isparallel) {
  library(parallel)
  # Determine optimal cluster size
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  
  # Export all variables to worker processes
  clusterExport(cl, ls())
  
  # Load required libraries on each worker process
  suppressMessages(clusterEvalQ(cl, expr = {
    library(ROI)     # R Optimization Infrastructure
    library(abess)   
    library(CVXR)    # Convex optimization (needed for LogRelax)
  }))
} else {
  # Load libraries for single-threaded execution
  library(ROI)
  library(abess)
}

# MAIN SIMULATION LOOP
# Iterate through each problem dimension
for (p in p_list) {
  # Track computation time for this dimension
  time_flag <- proc.time()[3][[1]]
  
  # RESULTS MATRIX INITIALIZATION
  # Create storage for results across all graph types
  # Note: 10 * res_num accounts for additional AUC metric (versus 8 in other simulations)
  res_summary <- matrix(0, nrow = length(type_list), ncol = 10 * res_num)
  colnames(res_summary) <- c(
    t(outer(method_name, rate_name, paste0)),    # Classification metrics
    t(outer(method_name, loss_name, paste0)),    # Reconstruction losses
    paste0(method_name, "prop"),                 # Exact recovery proportion
    paste0(method_name, "mse"),                  # Mean squared error
    paste0(method_name, "auc"),                  # Area under ROC curve
    "runtime"                                    # Computation time
  )
  
  # Storage for detailed results from each graph type
  res <- list()
  
  # GRAPH TYPE ITERATION
  # Test each graph structure type with current dimension
  for (k in 1:length(type_list)) {
    type <- type_list[k]
    print(paste0("p = ", p, "; type = ", type))
    
    # SIMULATION EXECUTION
    # Run multiple replications for statistical reliability
    if (isparallel) {
      # Parallel execution across worker processes
      clusterExport(cl, c("type", "n", "p", "method", "alpha", "beta", "degree"))
      res[[k]] <- parSapply(cl, 1:nrep, sim, type = type, n = n, p = p, method = method, 
                           alpha = alpha, beta = beta, degree = degree, thres = 0.0)
    } else {
      # Sequential execution (for debugging or small problems)
      res[[k]] <- sapply(1, sim, type = type, n = n, p = p, method = method, 
                        alpha = alpha, beta = beta, degree = degree, thres = 0.0)
    }
    
    # RESULTS AGGREGATION
    # Compute mean performance across all replications
    # Note: thres = 0.0 enables AUC computation (soft thresholding)
    res_summary[k, ] <- rowMeans(res[[k]], na.rm = TRUE)
  }
  
  # TIMING AND OUTPUT
  # Report computation time for this dimension
  print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))
  
  # SAVE RESULTS FOR THIS DIMENSION
  # Write comprehensive results to CSV file
  # Transpose for standard format (methods as rows, metrics as columns)
  write.csv(t(res_summary),
            paste0(path, result_file, method, "_p", p, "_", mark, ".csv"))
}

# CLEANUP
# Terminate parallel cluster if it was created
if (isparallel) {
  stopCluster(cl)
}