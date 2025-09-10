# ==============================================================================
# SIMULATION STUDY: WEAK SIGNAL REGIME
# ==============================================================================
# Purpose: This script evaluates graphical model estimation methods in the
# weak signal regime where interaction strengths are systematically varied
# from strong to weak. This tests method sensitivity and robustness when
# the signal-to-noise ratio decreases. The "ws" stands for "weak signal".
#
# Input: Command line arguments specifying method and graph type
# Output: CSV files containing minimum sample sizes across signal strength values
# ==============================================================================

# ENVIRONMENT INITIALIZATION
# Clean workspace and prepare environment
rm(list = ls()); gc(reset = TRUE)

# Set working directory and load simulation framework
path <- here::here(); setwd(path)
source("method_implementation.R")  # Graphical model estimation methods
source("evaluation.R")             # Performance evaluation functions
source("simulation_main.R")        # Main simulation control logic

# COMMAND LINE ARGUMENT PROCESSING
# Parse arguments to determine experimental configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Default settings when no arguments provided
  type_list <- c(1)      # Graph type: 1=random regular, 3=four-nearest neighbor
  method <- c("SLIDE")   # Estimation method to evaluate
} else {
  # Extract configuration from command line arguments
  for (arg in args) {
    if (grepl("^--method=", arg)) {
      method <- sub("^--method=", "", arg)
    }
    if (grepl("^--type=", arg)) {
      type_list <- as.integer(sub("^--type=", "", arg))
    }
  }
}

# SIMULATION PARAMETERS
# Core experimental settings
nrep <- 45              # Number of Monte Carlo replications
isparallel <- TRUE      # Enable parallel processing
ncore <- 45             # Number of CPU cores to utilize
save <- FALSE           # Whether to save intermediate results

# SAMPLE SIZE SEARCH CONFIGURATION
# Note: Sample sizes in code represent training set size
# For plotting, multiply by 2 to get total sample size
n_max <- 1e9            # Maximum sample size limit
p_list <- 16            # Fixed problem dimension

# GRAPH TYPE-SPECIFIC WEAK SIGNAL CONFIGURATION
# Different graph types require different signal strength ranges
if (type_list == 1) {
  ##### RANDOM REGULAR GRAPHS - WEAK SIGNAL REGIME #####
  n_start <- 200          # Starting sample size for search
  omega <- 1.2            # Total interaction strength constraint
  degree_list <- c(3)     # Fixed node degree
  
  # Create decreasing sequence of signal strengths (strong to weak)
  # alpha controls minimal signal strength: stronger alpha = easier recovery
  alpha_list <- 2 * seq(90, 9, length=10) / 450  # 10 points from strong to weak
  
  # Distribute remaining interaction strength across edges
  beta_list <- (omega - alpha_list) / (degree_list - 1)
  
} else if (type_list == 3) {
  ##### 4-NEAREST NEIGHBOR GRIDS - WEAK SIGNAL REGIME #####
  n_start <- 1000         # Higher starting point due to grid structure difficulty
  omega <- 0.9            # Lower total strength for grid graphs
  degree_list <- c(4)     # Fixed degree for 4-nearest neighbor
  
  # Adjusted signal strength range for grid graphs
  alpha_list <- 3 * (seq(30, 3, length=10) / 400)  # Different scaling for grids
  
  # Compute corresponding interaction strengths
  beta_list <- (omega - alpha_list) / (degree_list - 1)
}

# OUTPUT DIRECTORY
dir_name <- "result_ws"  # Directory for weak signal results
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
result_file <- paste0("/", dir_name, "/")

# PARAMETER COMBINATION SETUP
# Create matrix of all parameter combinations to test
l1 <- length(p_list); l2 <- length(alpha_list); 
l3 <- length(beta_list); l4 <- length(degree_list); 
conf_mat <- cbind(p_list, alpha_list, beta_list, degree_list)

# RESULTS ORGANIZATION
# Define structure for storing and naming results
method_name <- paste0(method, ".")           # Method name with separator
res_num <- length(method)                    # Number of methods (always 1)
rate_name <- c("TPR", "FPR", "MCC")         # Classification metrics
loss_name <- c("loss_op", "loss_l1", "loss_F")  # Reconstruction losses
mark <- "test"                               # Experiment identifier

# PARALLEL PROCESSING SETUP
# Initialize parallel computing if enabled
if (isparallel && nrep > 1) {
  library(parallel)
  # Use minimum of requested cores and replications
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  
  # Export all variables to worker processes
  clusterExport(cl, ls())
  
  # Load required libraries on each worker process
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
  
  # Iterate through each signal strength configuration
  for(conf_index in 1:nrow(conf_mat)) {
    conf <- conf_mat[conf_index, ]
    p <- conf[1]; alpha <- conf[2]; beta <- conf[3]; degree <- conf[4]
    
    # Create descriptive identifier for current configuration
    info <- paste0("p", p, "_alpha", alpha, "_beta", beta, "_degree", degree)
    print(paste0("config: ", info))
    
    # MINIMUM SAMPLE SIZE SEARCH FOR WEAK SIGNALS
    # Find smallest sample size achieving perfect recovery at this signal strength
    n_temp <- n_start_inner   # Current sample size being tested
    prop <- 0; k <- 1         # Recovery proportion and iteration counter
    
    # Initialize results storage matrix
    res_summary <- matrix(0, nrow = 1, ncol = 8 * res_num)
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), 
                               t(outer(method_name, loss_name, paste0)),
                               paste0(method_name, "prop"), 
                               paste0(method_name, "mse"))
    
    # Seed scaling for reproducible random number generation across iterations
    seed_scale <- 1
    
    # Sample size search loop - continue until perfect recovery or limit reached
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      
      # EXECUTE SIMULATION REPLICATIONS
      # Run multiple independent trials to estimate recovery probability
      if(isparallel == TRUE && nrep > 1) {
        # Parallel execution with scaled seeds for reproducibility
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, seed_scale * (1:nrep), sim, 
                          type = type, n = n_temp, p = p, method = method, 
                          alpha = alpha, beta = beta, degree = degree)
      } else {
        # Sequential execution (for debugging or single replication)
        temp <- sapply(1, sim, type = type, n = n_temp, p = p, method = method, 
                       alpha = alpha, beta = beta, degree = degree)
      }
      
      # RESULTS PROCESSING
      # Display individual recovery outcomes for monitoring weak signal performance
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
      # Optionally save progress for long-running weak signal experiments
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      # CHECK TERMINATION CONDITION
      # Stop search when perfect recovery achieved
      if(prop == 1) {
        # Set starting point for next signal strength slightly above current
        # This optimization works because weaker signals typically need larger samples
        n_start_inner <- round(n_temp * 1.005)
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      
      # ADAPTIVE SAMPLE SIZE UPDATE
      # Increase sample size based on current recovery rate
      # Formula: n_new = n_current * (2.0 - recovery_rate)
      # This provides finding the empirical sample complexity than linear search 
      # important for weak signals where sample size requirements can be very large
      n_temp <- round(n_temp * (2.0 - prop))
      k <- k + 1
    }
    
    # Handle computational limit exceeded
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary, rep(0, 8))  # Add failure indicator
    }
    
    print("______________________________________________________________")
    
    # SAVE FINAL RESULTS FOR THIS SIGNAL STRENGTH
    # Write complete results for this signal strength configuration
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
  }
}

# CLEANUP AND TIMING REPORT
# Shut down parallel cluster if it was created
if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

# Report total execution time for weak signal analysis
print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))