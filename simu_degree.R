# ==============================================================================
# SIMULATION STUDY: MAXIMUM NODE DEGREE VARIATION
# ==============================================================================
# Purpose: This script evaluates graphical model estimation methods across
# different maximum node degree values. Node degree affects the sparsity the graph, 
# which impacts the difficulty of structure recovery.
# Two scenarios are tested: local dense networks and random regular graphs.
#
# Input: Command line arguments specifying method and case type
# Output: CSV files containing minimum sample sizes for perfect recovery
# ==============================================================================

# ENVIRONMENT SETUP
# Clean workspace and prepare environment
rm(list = ls()); gc(reset = TRUE)

# Set working directory and load simulation framework
path <- here::here(); setwd(path)
source("method_implementation.R")  # Graphical model estimation methods
source("evaluation.R")             # Performance evaluation metrics
source("simulation_main.R")        # Main simulation control logic

# SIMULATION CONFIGURATION
# Core simulation parameters
nrep <- 45              # Number of Monte Carlo replications per configuration
isparallel <- TRUE      # Enable parallel processing for efficiency  
ncore <- 45             # Number of CPU cores to utilize
save <- FALSE           # Default: do not save intermediate results

# COMMAND LINE ARGUMENT PROCESSING
# Parse arguments to determine experimental scenario
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Default configuration when no arguments provided
  method <- c("SLIDE")   # Method to evaluate
  third_order <- FALSE   # Type of network structure
} else {
  # Extract method name from --method= argument
  method <- sub("^--method=", "", args[grep("^--method=", args)])
  # Extract case number from --case= argument to determine network type
  case <- as.numeric(sub("^--case=", "", args[grep("^--case=", args)]))
  if (case == 0) {
    third_order <- FALSE  # Random regular graphs
  } else {
    third_order <- TRUE   # Local dense networks
  }
}

# EXPERIMENTAL DESIGN PARAMETERS
# Maximum sample size for computational feasibility
n_max <- 1e9

# Fixed number of variables across all experiments
p_list <- c(16)

# SCENARIO-SPECIFIC PARAMETER CONFIGURATION
if (third_order) {
  ##### LOCAL DENSE NETWORK SCENARIO #####
  # Test networks with high local connectivity
  n_start <- 2e1          # Starting sample size (small due to easy recovery)
  type_list <- c(6)       # Graph type 6: local dense networks
  omega <- 2.0            # Total interaction strength constraint
  degree_list <- c(1:14)  # Range of node degrees to test
  
  # Parameter relationships for local dense networks
  # alpha + beta*degree = omega (interaction strength constraint)
  alpha_list <- omega / degree_list  # Minimal signal strength
  beta_list <- omega / degree_list   # Interaction strength
} else {
  ##### RANDOM REGULAR GRAPH SCENARIO #####
  # Test sparse random graphs with uniform degree
  n_start <- 5e4          # Starting sample size (larger due to harder recovery)
  omega <- 2.8            # Total interaction strength
  type_list <- c(7)       # Graph type 7: random regular graphs
  degree_list <- c(3:12)  # Degrees to test
  alpha_list <- c(0.2)    # Fixed Minimal signal strength
  
  # Distribute remaining interaction strength across edges
  # Formula: alpha + beta*(degree-1) = omega
  beta_list <- (omega - alpha_list) / (degree_list - 1)
}

# OUTPUT DIRECTORY SETUP
# Create directory
dir_name <- "result_d"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

# PARAMETER COMBINATION MATRIX
# Create all combinations of experimental parameters
l1 <- length(p_list); l2 <- length(alpha_list); 
l3 <- length(beta_list); l4 <- length(degree_list); 
conf_mat <- cbind(p_list, alpha_list, beta_list, degree_list)

# RESULTS STRUCTURE SETUP
# Define metric names and result organization
method_name <- paste0(method, ".")           # Add separator to method name
res_num <- length(method)                    # Number of methods being tested
rate_name <- c("TPR", "FPR", "MCC")         # Classification accuracy metrics
loss_name <- c("loss_op", "loss_l1", "loss_F")  # Matrix reconstruction losses
mark <- "test"                               # Experiment identifier

# PARALLEL PROCESSING INITIALIZATION
# Set up parallel computing cluster if requested
if (isparallel && nrep > 1) {
  library(parallel)
  # Determine optimal number of worker processes
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  
  # Share all variables with worker processes
  clusterExport(cl, ls())
  
  # Load required packages on each worker
  suppressMessages(clusterEvalQ(cl, expr = {
    library(ROI)     
    library(abess)    
  }))
}

# MAIN SIMULATION EXECUTION
# Track total computation time for performance monitoring
time_flag <- proc.time()[3][[1]]

# Process each graph type
for(type in type_list) {
  print(paste0("type: ", type))
  n_start_inner <- n_start  # Initialize sample size search starting point
  
  # Test each parameter configuration
  for(conf_index in 1:nrow(conf_mat)) {
    conf <- conf_mat[conf_index, ]
    p <- conf[1]; alpha <- conf[2]; beta <- conf[3]; degree <- conf[4]
    
    # Generate descriptive name for current configuration
    info <- paste0("p", p, "_alpha", alpha, "_beta", beta, "_degree", degree)
    print(paste0("config: ", info))
    
    # MINIMUM SAMPLE SIZE SEARCH
    # Adaptively find smallest sample size achieving perfect recovery
    n_temp <- n_start_inner   # Current sample size being tested
    prop <- 0; k <- 1         # Recovery proportion and iteration number
    
    # Initialize results storage matrix
    res_summary <- matrix(0, nrow = 1, ncol = 8 * res_num)
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), 
                               t(outer(method_name, loss_name, paste0)), 
                               paste0(method_name, "prop"), 
                               paste0(method_name, "mse"))
    
    # Seed scaling for reproducible random number generation across iterations
    seed_scale <- 1
    
    # Sample size search loop
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      
      # EXECUTE SIMULATION REPLICATIONS
      # Run multiple replications to estimate recovery probability
      if(isparallel == TRUE && nrep > 1) {
        # Parallel execution with scaled seeds for reproducibility
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, seed_scale * (1:nrep), sim, type = type, n = n_temp, p = p, method = method, 
                          alpha = alpha, beta = beta, degree = degree)
      } else {
        # Single-threaded execution
        temp <- sapply(1, sim, type = type, n = n_temp, p = p, method = method, 
                       alpha = alpha, beta = beta, degree = degree)
      }
      
      # RESULTS PROCESSING
      # Display individual recovery outcomes for monitoring
      print(paste0("is recovery: ", paste0(temp[7, ], collapse = " ")))
      
      # Compute mean performance across all replications
      temp <- rowMeans(temp, na.rm = TRUE)
      prop <- temp[7]  # Extract exact recovery proportion
      
      # Store results and update tracking variables
      if(k == 1) {
        res_summary[k, ] <- temp
        seed_scale <- 1    # Reset seed scaling for first iteration
      } else {
        res_summary <- rbind(res_summary, temp)
        seed_scale <- seed_scale + 1  # Increment for different random streams
      }
      rownames(res_summary)[k] <- paste0("n_", n_temp)
      
      # SAVE INTERMEDIATE RESULTS
      # Optionally save progress for long-running experiments
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      # CHECK TERMINATION CONDITION
      # Stop search when perfect recovery achieved
      if(prop == 1) {
        # Set starting point for next configuration slightly above current
        n_start_inner <- round(n_temp * 1.005)
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      
      # ADAPTIVE SAMPLE SIZE UPDATE
      # Increase sample size based on current recovery rate
      n_temp <- round(n_temp * (2.0 - prop))
      k <- k + 1
    }
    
    # Handle computational limit exceeded
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary, rep(0, 8))  # Add failure indicator
    }
    
    print("______________________________________________________________")
    
    # SAVE FINAL RESULTS
    # Write complete results for this parameter configuration
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
  }
}

# CLEANUP AND TIMING REPORT
# Shut down parallel cluster if it was created
if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

# Report total execution time for performance analysis
print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))