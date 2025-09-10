# ==============================================================================
# SIMULATION STUDY: PROBLEM DIMENSION SCALING
# ==============================================================================
# Purpose: This script evaluates how graphical model estimation methods scale
# with increasing problem dimension (number of variables p). It finds the minimum
# sample size required for perfect recovery as p increases, testing scalability
# of different estimation methods across two graph types: 4-nearest neighbor
# grids and random regular graphs.
#
# Input: Command line arguments specifying graph type and estimation method
# Output: CSV files containing minimum sample sizes for each dimension
# ==============================================================================

# ENVIRONMENT SETUP
# Initialize clean workspace
rm(list = ls()); gc(reset = TRUE)

# Set working directory and load simulation framework
path <- here::here(); setwd(path)
source("method_implementation.R")  # Graphical model estimation methods
source("evaluation.R")             # Performance evaluation functions
source("simulation_main.R")        # Core simulation framework

# COMMAND LINE ARGUMENT PROCESSING
# Parse arguments to determine experimental configuration
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Default settings when no arguments provided
  type_list <- c(1)      # Graph type: 1=random regular, 3=4-nearest neighbor
  method <- c("SLIDE")   # Estimation method to evaluate
} else {
  # Extract configuration from command line arguments
  for (arg in args) {
    if (grepl("^--type=", arg)) {
      type_list <- as.integer(sub("^--type=", "", arg))
    }
    if (grepl("^--method=", arg)) {
      method <- sub("^--method=", "", arg)
    }
  }
}

# SIMULATION PARAMETERS
# Core experimental settings
nrep <- 45              # Number of Monte Carlo replications per configuration
isparallel <- TRUE      # Enable parallel computing for efficiency
ncore <- 45             # Number of CPU cores to use
save <- FALSE           # Whether to save intermediate results

# SAMPLE SIZE SEARCH BOUNDS
# Maximum sample size to prevent infinite search
n_max <- 1e9

# GRAPH TYPE-SPECIFIC CONFIGURATION
# Different graph types require different parameter ranges and starting points
if (type_list %in% c(3)) {
  ##### 4-NEAREST NEIGHBOR GRID GRAPHS #####
  # Grid graphs: each node connected to 4 neighbors in 2D lattice
  p_list <- c(9, 16, 25, 36, 49, 64, 81, 100)  # Perfect squares for grid structure
  # Method-specific starting sample sizes to speedup
  if (method == 'SLIDE') {
    n_start <- 1674      
  } else if (method %in% c('RPLE_thres', 'RISE_thres', "logRISE_thres")) {
    n_start <- 800
  } else if (method == "ELASSO_thres") {
    n_start <- 1000
  }
  degree_list <- c(4)    # Fixed degree for 4-nearest neighbor
  
} else if (type_list %in% c(1)) {
  ##### RANDOM REGULAR GRAPHS #####
  p_list <- c(8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52)
  # Method-specific starting sample sizes  to speedup
  if (method %in% c('SLIDE', 'RPLE_thres', 'RISE_thres', "logRISE_thres")) {
    n_start <- 800
  } else if (method == "ELASSO_thres") {
    n_start <- 1800
  }
  degree_list <- c(3)    # Fixed degree for random regular graphs
}

# INTERACTION STRENGTH CONFIGURATION
# Set up parameter relationships for consistent difficulty across dimensions
omega <- 1.5            # Total interaction strength constraint
alpha_list <- 0.3       # Minimal signal strength

# Distribute remaining interaction strength across edges
# Formula: alpha + beta*(degree-1) = omega  
beta_list <- (omega - alpha_list) / (degree_list - 1)

# OUTPUT DIRECTORY
dir_name <- "result_p"  # Directory for dimension scaling results
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
res_num <- length(method)                    # Number of methods (always 1 here)
rate_name <- c("TPR", "FPR", "MCC")         # Classification performance metrics
loss_name <- c("loss_op", "loss_l1", "loss_F")  # Matrix reconstruction losses
mark <- "test"                               # Experiment identifier

# PARALLEL PROCESSING SETUP
# Initialize parallel computing cluster if enabled
if (isparallel && nrep > 1) {
  library(parallel)
  # Use minimum of requested cores and replications
  cl.cores <- min(ncore, nrep)
  cl <- makeCluster(cl.cores)
  
  # Export all variables to worker processes
  clusterExport(cl, ls())
  
  # Load required libraries on each worker
  suppressMessages(clusterEvalQ(cl, expr = {
    library(ROI)     
    library(abess)   
  }))
}

# MAIN SIMULATION EXECUTION
# Track total computation time
time_flag <- proc.time()[3][[1]]

# Process each graph type (usually just one)
for(type in type_list) {
  print(paste0("type: ", type))
  n_start_inner <- n_start  # Initialize search starting point
  
  # Iterate through each problem dimension
  for(conf_index in 1:nrow(conf_mat)) {
    conf <- conf_mat[conf_index, ]
    p <- conf[1]; alpha <- conf[2]; beta <- conf[3]; degree <- conf[4]
    
    # Create descriptive identifier for current configuration
    info <- paste0("p", p, "_alpha", alpha, "_beta", beta, "_degree", degree)
    print(paste0("config: ", info))
    
    # MINIMUM SAMPLE SIZE SEARCH
    # Adaptive search to find smallest n achieving perfect recovery
    n_temp <- n_start_inner   # Current sample size being tested
    prop <- 0; k <- 1         # Recovery proportion and iteration counter
    
    # Initialize results storage
    res_summary <- matrix(0, nrow = 1, ncol = 8 * res_num)
    colnames(res_summary) <- c(t(outer(method_name, rate_name, paste0)), 
                               t(outer(method_name, loss_name, paste0)),
                               paste0(method_name, "prop"), 
                               paste0(method_name, "mse"))
    
    # Seed scaling for reproducible random number generation
    seed_scale <- 1
    
    # Sample size search loop
    while(n_temp <= n_max) {
      print(paste0("n_temp: ", n_temp))
      
      # EXECUTE SIMULATION REPLICATIONS
      # Run multiple independent trials to estimate recovery probability
      if(isparallel == TRUE && nrep > 1) {
        # Parallel execution with scaled seeds
        clusterExport(cl, c("type", "n_temp", "p", "method", "alpha", "beta", "degree"))
        temp <- parSapply(cl, seed_scale * (1:nrep), sim, 
                          type = type, n = n_temp, p = p, method = method, 
                          alpha = alpha, beta = beta, degree = degree)
      } else {
        # Sequential execution
        temp <- sapply(1, sim, type = type, n = n_temp, p = p, method = method, 
                       alpha = alpha, beta = beta, degree = degree)
      }
      
      # PROCESS SIMULATION RESULTS
      # Display recovery outcomes for monitoring
      print(paste0("is recovery: ", paste0(temp[7, ], collapse = " ")))
      
      # Compute average performance across replications
      temp <- rowMeans(temp, na.rm = TRUE)
      prop <- temp[7]  # Extract exact recovery proportion
      
      # Store results in summary matrix
      if(k == 1) {
        res_summary[k, ] <- temp
        seed_scale <- 1      # Reset seed scaling
      } else {
        res_summary <- rbind(res_summary, temp)
        seed_scale <- seed_scale + 1  # Increment for different random streams
      }
      rownames(res_summary)[k] <- paste0("n_", n_temp)
      
      # SAVE INTERMEDIATE RESULTS
      # Optional saving for long-running experiments
      if(save) {
        write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_",
                                         info, "_n", n_temp, "_", mark,".csv"))
      }
      
      print(paste0("Sample size = ", n_temp, ";  Recovery prop = ", prop))
      
      # CHECK TERMINATION CONDITION
      # Stop when perfect recovery achieved
      if(prop == 1) {
        n_start_inner <- n_temp + 1  # Set starting point for next dimension
        print(paste0("Sufficient size for config ", info, " = ", n_temp))
        break
      } 
      
      # ADAPTIVE SAMPLE SIZE UPDATE
      # Increase sample size using dimension-specific strategies
      n_increase_scale <- 1
      if (method == "ELASSO_thres") {
        n_increase_scale <- 4    # ELASSO needs larger increases
      }
      
      # Graph type-specific increment strategies
      if (type == 3) {
        # 4-nearest neighbor: use adaptive increment based on current sample size
        n_temp <- n_temp + max(round(n_increase_scale * 65536 / (n_temp + 1)), 1)
      } else if (type == 1) {
        # Random regular: smaller base increment
        n_temp <- n_temp + max(round(n_increase_scale * 16384 / (n_temp + 1)), 1)
      }
      
      k <- k + 1
    }
    
    # Handle computational limit exceeded
    if(n_temp > n_max) {
      print(paste0("Sufficient size for config ", info, " = more than n_max!"))
      res_summary <- rbind(res_summary, rep(0, 8))  # Add failure indicator
    }
    
    print("______________________________________________________________")
    
    # SAVE FINAL RESULTS
    # Write complete results for this dimension
    write.csv(t(res_summary), paste0(path, result_file, method, "_type", type, "_", 
                                     info, "_", mark,".csv"))
  }
}

# CLEANUP AND TIMING REPORT
# Shut down parallel cluster if created
if(isparallel == TRUE && nrep > 1) {
  stopCluster(cl)
}

# Report total execution time
print(paste0("time:", proc.time()[3][[1]] - time_flag, ' s'))