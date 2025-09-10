# ==============================================================================
# SIMULATION STUDY: MAXIMUM NEIGHBORHOOD WEIGHT
# ==============================================================================
# Purpose: This script conducts simulation studies to evaluate graphical model
# estimation methods across different max neighborhood weight values. The 
# parameter is controlled by the strength of (non-minimal) interactions (denoted  
# as beta) in the graphical model. The simulation finds the minimum sample size   
# required for perfect recovery for each beta value.
#
# Input: Command line arguments specifying method and graph type
# Output: CSV files containing performance metrics and minimum sample sizes
# ==============================================================================

# ENVIRONMENT SETUP
# Clear workspace and reset garbage collection for clean start
rm(list = ls()); gc(reset = TRUE)

# Load required libraries
library(stringr)  # String manipulation functions

# Set working directory and source required functions
path <- here::here(); setwd(path)
source("method_implementation.R")  # Graphical model estimation methods
source("evaluation.R")             # Performance evaluation functions  
source("simulation_main.R")        # Main simulation framework

# COMMAND LINE ARGUMENT PARSING
# Parse command line arguments to determine which method and graph type to test
# Arguments format: --method=METHOD_NAME --type=TYPE_NUMBER
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Default settings when no arguments provided
  # method can be: "RPLE_thres", "RISE_thres", "logRISE_thres", "ELASSO_thres", "SLIDE"
  method <- c("SLIDE")      # Method to evaluate
  type_list <- c(1)         # Graph type to test
} else {
  # Extract method name from --method= argument
  method <- sub("^--method=", "", args[grep("^--method=", args)])
  # Extract graph type from --type= argument  
  type_list <- as.numeric(sub("^--type=", "", args[grep("^--type=", args)]))
}

# SIMULATION PARAMETERS
# Core simulation settings
nrep <- 45              # Number of repetitions for each configuration
isparallel <- TRUE      # Enable parallel computing for efficiency
ncore <- 45             # Number of CPU cores to use
save <- FALSE           # Whether to save intermediate results

# EXPERIMENTAL DESIGN PARAMETERS
# Fixed parameters across all simulations
p_list <- c(16)         # Number of variables (nodes) in the graph
alpha_list <- c(0.4)    # Minimal signal strength parameter
degree_list <- c(3)     # Average node degree in the graph

# SAMPLE SIZE SEARCH PARAMETERS
# Adaptive sample size search to find minimum required for perfect recovery
n_start <- 45711          # Starting sample size for search
if (str_detect(method, "LogRelax")) {
  n_max <- 1e6          # Maximum sample size for LogRelax (computationally intensive)
} else {
  n_max <- 2e9          # Maximum sample size for other methods
}

# BETA PARAMETER CONFIGURATION
# Beta controls higher-order interaction strength
# Different graph types require different beta ranges for meaningful variation
if (type_list == 1) {
  # Type B in the main paper
  beta_list <- (14:24)
} else if (type_list == 2) {
  # Type D in the main paper
  beta_list <- (20:32)
} else if (type_list == 3) {
  # Type A in the main paper
  beta_list <- (20:36) / 2
} else if (type_list == 4) {
  # Type C in the main paper
  beta_list <- (22:32)
} else if (type_list == 5) {
  # Type E in the main paper
  beta_list <- (20:36) / 2
}
beta_list <- beta_list / 20

# OUTPUT DIRECTORY SETUP
# Create type-specific directory for results
dir_name <- paste0("result_type", type_list)
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

# PARAMETER COMBINATION SETUP
# Sort parameter lists for systematic exploration
p_list <- sort(p_list)
alpha_list <- sort(alpha_list, decreasing = TRUE)
beta_list <- sort(beta_list)
degree_list <- sort(degree_list)

# Create configuration matrix with all parameter combinations
l1 <- length(p_list); l2 <- length(alpha_list); l3 <- length(beta_list); l4 <- length(degree_list); 
conf_mat <- cbind(p_list, alpha_list, beta_list, degree_list)

# RESULTS STRUCTURE SETUP
# Define naming conventions for results
method_name <- paste0(method, ".")   # Method name with separator
res_num <- length(method)            # Number of methods being compared
rate_name <- c("TPR", "FPR", "MCC")  # Classification performance metrics
loss_name <- c("loss_op", "loss_l1", "loss_F")  # Matrix reconstruction losses
mark <- "test"                       # Identifier for this simulation run

# PARALLEL COMPUTING SETUP
# Configure parallel processing if enabled and multiple replications requested
if (isparallel && nrep > 1) {
  library(parallel)
  # Use minimum of requested cores and available replications
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  
  # Export all variables to worker processes
  clusterExport(cl, ls())
  
  # Load required libraries on each worker process
  suppressMessages(clusterEvalQ(cl, expr = {
    library(ROI)
    library(abess)
    library(CVXR)
  }))
}

# MAIN SIMULATION LOOP
# Track total computation time
time_flag <- proc.time()[3][[1]]
stop_flag <- FALSE

# Iterate through each graph type
for(type in type_list) {
  print(paste0("type: ", type))
  n_start_inner <- n_start  # Initialize sample size for this type
  
  # Iterate through each parameter configuration
  for(conf_index in 1:nrow(conf_mat)) {
    conf <- conf_mat[conf_index, ]
    p <- conf[1]; alpha <- conf[2]; beta <- conf[3]; degree <- conf[4]
    
    # Create descriptive identifier for this configuration
    info <- paste0("p", p, "_alpha", alpha, "_beta", beta, "_degree", degree)
    print(paste0("config: ", info))
    
    # ADAPTIVE SAMPLE SIZE SEARCH
    # Find minimum sample size required for perfect recovery
    n_temp <- n_start_inner   # Starting sample size for this configuration
    prop <- 0; k <- 1         # Recovery proportion and iteration counter
    
    # Initialize results matrix for this configuration
    res_summary <- matrix(0, nrow = 1, ncol = 8 * res_num)
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), 
                               t(outer(method_name, loss_name, paste0)),
                               paste0(method_name, "prop"), 
                               paste0(method_name, "mse"))
    
    # Sample size search loop - continue until perfect recovery or max sample size
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      
      # RUN SIMULATION REPLICATIONS
      # Execute simulation with current sample size
      if(isparallel == TRUE && nrep > 1) {
        # Parallel execution across multiple cores
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, 1:nrep, sim, type = type, n = n_temp, p = p, method = method, 
                          alpha = alpha, beta = beta, degree = degree)
      } else {
        # Sequential execution (single core or single replication)
        temp <- sapply(1, sim, type = type, n = n_temp, p = p, method = method, 
                       alpha = alpha, beta = beta, degree = degree)
      }
      
      # PROCESS SIMULATION RESULTS
      # Display recovery results for diagnostics
      print(paste0("is recovery: ", paste0(temp[7, ], collapse = " ")))
      
      # Compute average performance across replications
      temp <- rowMeans(temp, na.rm = TRUE)
      prop <- temp[7]  # Exact recovery proportion
      
      # Store results in summary matrix
      if(k == 1) {
        res_summary[k, ] <- temp
      } else {
        res_summary <- rbind(res_summary, temp)
      }
      rownames(res_summary)[k] <- paste0("n_", n_temp)
      
      # SAVE RESULTS
      # Save current results if requested
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      # CHECK TERMINATION CONDITIONS
      # Stop if perfect recovery achieved
      if(prop == 1) {
        n_start_inner <- n_temp  # Update starting point for next configuration
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      
      # Increase sample size based on current recovery rate
      # Formula: n_new = n_current * (2 - recovery_rate)
      # This provides faster finding the empirical sample complexity than linear search
      n_temp <- round(n_temp * (2 - prop))
      k <- k + 1
    }
    
    # Handle case where maximum sample size exceeded
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary,rep(0, 8))  # Add zero row to indicate failure
      stop_flag <- TRUE
    }
    
    print("______________________________________________________________")
    
    # SAVE FINAL RESULTS FOR THIS CONFIGURATION
    # Write complete results for this parameter configuration
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
    
    # Stop if maximum sample size exceeded (computational limit reached)
    if (stop_flag) {
      print("Stop sufficient sample size searching!")
      break;
    }
  }
}

# CLEANUP AND FINAL REPORTING
# Stop parallel cluster if it was created
if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

# Report total computation time
print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))