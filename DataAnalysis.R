# ==============================================================================
# US SENATE VOTING NETWORK ANALYSIS
# ==============================================================================
# Purpose: This script analyzes US Senate voting patterns using graphical model
# estimation to identify political networks and coalitions. It processes voting
# records from multiple Congressional sessions (112th-117th Congress, 2011-2023)
# and applies the SLIDE method to estimate sparse precision matrices representing
# senator relationships based on voting similarity.
#
# Input: 
#   - Sall_members.csv: Senator demographic and party information
#   - Sall_votes.csv: Individual voting records by senator and roll call
#
# Output:
#   - Network visualizations showing political coalitions
#   - Saved data objects for further analysis
# ==============================================================================

# ENVIRONMENT SETUP
# Clear workspace and prepare clean environment for analysis
rm(list = ls()); gc(reset = TRUE)

# Set working directory using here package for reproducible file paths
path <- here::here(); setwd(path)

# LOAD REQUIRED LIBRARIES
# Data manipulation and visualization packages
library(reshape2)    # Data reshaping functions (dcast, melt)
library(dplyr)       # Data manipulation and pipeline operations
library(abess)       # The SLIDE estimator
library(ggnetwork)   # Network visualization with ggplot2
library(ggpmisc)     # Additional ggplot2 extensions
library(network)     # Network object creation and manipulation
library(sna)         # Social network analysis functions
library(ggpubr)      # Publication-ready ggplot2 themes

# STUDY PERIOD DEFINITION
# Define the Congressional sessions and years to analyze
congress_study <- 112:117  # 112th through 117th Congress
year <- seq(2011, 2023, by = 2)  # Corresponding years (each Congress spans 2 years)

# ==============================================================================
# SENATE MEMBER INFORMATION PROCESSING
# ==============================================================================
# Load and process senator demographic and party affiliation data

# LOAD SENATOR INFORMATION
# Read CSV file containing senator biographical and political information
senate_info <- read.csv("Sall_members.csv")

# FILTER FOR STUDY PERIOD
# Keep only senators from the Congressional sessions of interest
senate_info <- senate_info[senate_info[["congress"]] %in% congress_study, ]

# STANDARDIZE PARTY CODING
# Convert party codes to standard format: 100=Democrat, 200=Republican, others=Independent
senate_info[["party_code"]][!(senate_info[["party_code"]] %in% c(100, 200))] <- 300

# CREATE READABLE PARTY LABELS
# Convert numeric party codes to letter abbreviations for easier interpretation
senate_info[["party"]][senate_info[["party_code"]] == 100] <- "D"  # Democrat
senate_info[["party"]][senate_info[["party_code"]] == 200] <- "R"  # Republican
senate_info[["party"]][senate_info[["party_code"]] == 300] <- "I"  # Independent

# COMPUTE DERIVED VARIABLES
# Calculate the calendar year for each Congress (starts in odd years)
senate_info[["year"]] <- 2 * senate_info[["congress"]] + 1789

# Calculate senator age during their service
senate_info[["age"]] <- senate_info[["year"]] - senate_info[["born"]]

# VERIFY PARTY DISTRIBUTION
# Display party composition across all studied Congressional sessions
table(senate_info[["party"]])

# SELECT RELEVANT COLUMNS AND SPLIT BY CONGRESS
# Keep only essential variables for network analysis
senate_info <- senate_info[c("congress", "chamber", "icpsr", "bioname", "party", "age")]

# Split data by Congressional session for separate analysis of each period
# This creates a list where each element contains senators for one Congress
senate_info <- split(senate_info, f = senate_info[["congress"]])

# ==============================================================================
# VOTING RECORDS PROCESSING
# ==============================================================================
# Load and process individual voting records to create senator-vote matrices

# LOAD VOTING DATA
# Read CSV file containing individual votes by senator and roll call number
votes <- read.csv("Sall_votes.csv")

# FILTER FOR STUDY PERIOD
# Keep only votes from the Congressional sessions of interest
votes <- votes[votes[["congress"]] %in% congress_study, ]

# VERIFY DATA STRUCTURE
# Check that we only have Senate votes (not House)
unique(votes[["chamber"]])

# REMOVE UNNECESSARY COLUMNS
# Remove chamber and probability columns not needed for network analysis
votes[["chamber"]] <- NULL
votes[["prob"]] <- NULL

# EXAMINE AND STANDARDIZE VOTE CODING
# Display distribution of different vote types
table(votes[["cast_code"]])

# BINARY VOTE ENCODING
# Convert vote codes to binary: 1 for "Yes" votes, -1 for all others (No/Abstain)
# This simplifies the analysis to focus on agreement patterns
votes[["cast_code"]] <- ifelse(votes[["cast_code"]] == 1, 1, -1)

# SPLIT BY CONGRESSIONAL SESSION
# Create separate voting matrices for each Congress
votes <- split(votes, f = votes[["congress"]])

# CREATE SENATOR-VOTE MATRICES
# Transform long-format voting data into wide-format matrices
votes <- lapply(votes, function(x) {
  # RESHAPE DATA
  # Create matrix: rows=senators (icpsr), columns=roll call votes, values=vote choices
  vote_x <- dcast(x, formula = icpsr ~ rollnumber, value.var = "cast_code")
  
  # SET ROW NAMES
  # Use senator ICPSR codes as row identifiers
  rownames(vote_x) <- vote_x[["icpsr"]]
  vote_x[["icpsr"]] <- NULL
  
  # REMOVE SENATORS WITH EXCESSIVE MISSING VOTES
  # Exclude senators who missed more than 0% of votes (currently keeps all)
  exclude_index <- which(rowMeans(is.na(vote_x)) > 0.0)
  if (length(exclude_index) != 0) {
    vote_x <- vote_x[-exclude_index, ]
  }
  
  # REMOVE UNANIMOUS VOTES
  # Exclude roll calls where all senators voted the same way
  # These provide no information for distinguishing voting patterns
  exclude_index <- which(apply(vote_x, 2, function(x) {
    length(unique(x))
  }) == 1)
  if (length(exclude_index) != 0) {
    vote_x <- vote_x[, -exclude_index]
  }
  
  # Return cleaned voting matrix
  vote_x
})

# VERIFY DATA QUALITY
# Check for remaining missing values and matrix dimensions
sapply(votes, anyNA)    # Should be FALSE for all Congressional sessions
sapply(votes, dim)      # Display number of senators and votes per Congress

# ==============================================================================
# DATA INTEGRATION AND PREPARATION
# ==============================================================================
# Combine senator information with voting matrices for network analysis

# MERGE SENATOR INFO WITH VOTING DATA
# Create integrated dataset combining demographics with voting matrices
senate_model_data <- mapply(function(x, y) {
  # ALIGN SENATOR INFORMATION WITH VOTING MATRIX
  # Match senators in voting matrix with their demographic information
  index <- match(as.numeric(rownames(y)), x[["icpsr"]])
  x <- x[index, ]
  
  # CLEAN UP VARIABLE NAMES
  # Remove ICPSR identifier column and rename biographical name field
  x[["icpsr"]] <- NULL
  colnames(x)[colnames(x) == "bioname"] <- "name"
  
  # RETURN STRUCTURED DATA
  # Create list with senator information and voting matrix as separate components
  list("info" = x, "vote" = as.matrix(y))
}, senate_info, votes, SIMPLIFY = FALSE)

# SAVE INTEGRATED DATA
# Store processed data for reproducible analysis
save(senate_model_data, file = "senate_model_data.rda")

# DISPLAY DATA SUMMARY
# Show dimensions of voting data for each Congressional session
tab_info <- sapply(senate_model_data, function(x) {
  dim(t(x[["vote"]]))  # Transpose to show votes x senators
})
print(tab_info)

# ==============================================================================
# NETWORK ESTIMATION USING SLIDE METHOD
# ==============================================================================
# Apply SLIDE (Sparse Learning with Interaction Detection and Estimation)
# to identify political networks based on voting similarity patterns

# LOAD PROCESSED DATA
load("senate_model_data.rda")

# ESTIMATE SPARSE PRECISION MATRICES
# Apply SLIDE method to each Congressional session separately
senate_nodewise <- lapply(senate_model_data, function(x) {
  # PREPARE DATA FOR SLIDE
  # Transpose voting matrix so variables (senators) are columns
  dat <- t(x[["vote"]])
  
  # APPLY SLIDE ALGORITHM
  # Parameters:
  #   - dat: voting data matrix (votes x senators)
  #   - weight: equal weighting for all votes
  #   - tune.type: use Generalized Information Criterion for model selection
  #   - ic.scale: penalty scaling factor (1 = standard BIC penalty)
  splicing_gic_logistic <- slide(dat, weight = rep(1, nrow(dat)), 
                                 tune.type = "gic", ic.scale = 1)
  
  # EXTRACT ESTIMATED PRECISION MATRIX
  # Return the estimated sparse precision matrix representing senator relationships
  splicing_gic_logistic[[1]]
})

# SAVE NETWORK ESTIMATES
# Store estimated precision matrices for visualization and analysis
save(senate_nodewise, file = "senate_nodewise.rda")

# ==============================================================================
# NETWORK VISUALIZATION FUNCTION
# ==============================================================================
# Function to convert sparse precision matrices into network plot data

#' Convert Precision Matrix to Network Plot Data
#' 
#' Transforms a sparse precision matrix into a format suitable for network
#' visualization using ggnetwork and ggplot2.
#' 
#' @param theta Sparse precision matrix representing senator relationships
#' @param dat Data frame with senator information (name, party, etc.)
#' @param complete Logical; if TRUE, include all nodes even if isolated
#' @return ggnetwork object ready for plotting with ggplot2
#' 
#' Output contains:
#   - x, y: node coordinates for plotting
#   - xend, yend: edge endpoint coordinates
#   - connection_strength: edge weights based on precision matrix values
#   - color: node colors based on party affiliation
#   - vertex.names: senator names for labeling
plot_data <- function(theta, dat, complete = FALSE) {
  # IDENTIFY ACTIVE EDGES
  # Find non-zero elements in precision matrix (significant relationships)
  act_set <- which(theta != 0, arr.ind = TRUE)
  
  # KEEP ONLY UPPER TRIANGLE
  # Avoid duplicate edges since precision matrix is symmetric
  act_set <- act_set[act_set[, 1] < act_set[, 2],]
  
  # CONVERT INDICES TO SENATOR NAMES
  # Replace matrix indices with senator names for interpretability
  tmp <- act_set
  for(i in 1:nrow(act_set)) {
    for(j in 1:ncol(act_set)) {
      act_set[i, j] <- colnames(theta)[tmp[i, j]]
    }
  }
  act_set <- apply(act_set, 2, as.character)
  
  # CALCULATE NETWORK PROPERTIES
  p <- ncol(theta)                                      # Total number of senators
  p1 <- length(unique(as.vector(act_set)))             # Number of connected senators
  congress <- dat[["congress"]][1]                      # Congressional session number
  nodesize <- 4                                        # Default node size for plotting
  
  # CREATE NETWORK OBJECT
  if(complete) {
    # CREATE COMPLETE NETWORK (including isolated nodes)
    g <- network::network.initialize(p, directed = FALSE)
    network::network.edgelist(tmp, g)
    network::network.vertex.names(g) <- colnames(theta)
  } else {
    # CREATE NETWORK WITH ONLY CONNECTED NODES
    g <- network::network(act_set, matrix.type = "edgelist", directed = FALSE)
  }
  
  # CALCULATE EDGE WEIGHTS
  # Use absolute values of precision matrix elements as connection strengths
  connection_strength <- vector(length = nrow(act_set))
  for(i in 1:length(connection_strength)) {
    connection_strength[i] <- abs(theta[tmp[i, 1], tmp[i, 2]])
  }
  connection_strength <- as.numeric(connection_strength)
  
  # ASSIGN EDGE ATTRIBUTES
  # Store connection strength as network edge attribute
  g %v% "connection_strength" <- connection_strength
  
  # ASSIGN NODE ATTRIBUTES
  # Match network vertices with senator information
  vertex_name <- network::network.vertex.names(g)
  match_id <- match(vertex_name, dat$name)
  party <- dat$party[match_id]
  
  # Set node colors based on party affiliation
  g %v% "color" <- as.character(party)
  g %v% "text" <- ""  # Initialize empty text labels
  
  # CONVERT TO GGNETWORK FORMAT
  # Transform network object to data frame format for ggplot2
  # Uses eigenvector-based layout for node positioning
  ggnetwork_default <- ggnetwork::ggnetwork(g, layout = "eigen")
  
  return(ggnetwork_default)
}

# ==============================================================================
# VISUALIZATION DATA PREPARATION
# ==============================================================================
# Load network estimates and prepare data for multi-panel visualization

# LOAD SAVED NETWORK ESTIMATES AND ORIGINAL DATA
load("senate_nodewise.rda")      # Estimated precision matrices
load("senate_model_data.rda")    # Original senator and voting data

# FILTER DATA FOR VISUALIZATION
# Select specific Congressional sessions for plotting
congress_range <- 112:117
index <- which(names(senate_nodewise) %in% as.character(congress_range))
senate_nodewise <- senate_nodewise[index]
senate_model_data <- senate_model_data[index]

# CREATE PLOT DATA FOR ALL CONGRESSIONAL SESSIONS
# Generate network plot data for each Congress with appropriate labeling
pdat_list <- mapply(function(x, y) {
  # PREPARE PRECISION MATRIX
  # Extract senator information and use names as matrix column names
  y <- y[["info"]]
  colnames(x) <- y[["name"]]
  
  # CONVERT TO PLOT DATA
  # Transform precision matrix to ggnetwork format
  pdat <- plot_data(x, y, complete = FALSE)
  
  # ADD CONGRESSIONAL SESSION LABELS
  # Calculate start and end years for this Congress
  start_year <- 2 * y[["congress"]][1] + 1787
  end_year <- start_year + 2
  congress_info <- paste0(sprintf("%s-th Congress", y[["congress"]][1]), 
                          " (", start_year, "-", end_year, ")")
  
  # ADD CONGRESS IDENTIFIER TO PLOT DATA
  pdat <- cbind(pdat, "congress" = congress_info)
  pdat
}, senate_nodewise, senate_model_data, SIMPLIFY = FALSE)

# COMBINE ALL CONGRESSIONAL SESSIONS
# Merge plot data from all sessions into single data frame for faceted plotting
pdat <- do.call("rbind.data.frame", pdat_list)

# ==============================================================================
# NETWORK VISUALIZATION
# ==============================================================================
# Create publication-quality network plots showing political coalitions

# SET VISUALIZATION PARAMETERS
nodesize <- 3        # Node size for plotting
set.seed(123)       # Set seed for reproducible node layouts

# CREATE BASIC NETWORK PLOT
# Multi-panel visualization showing political networks across Congressional sessions
basic_pic <- ggplot2::ggplot(pdat, aes(x, y, xend = xend, yend = yend)) +
  facet_wrap(congress ~ ., scales = "free", ncol = 3) + # CREATE SEPARATE PANELS FOR EACH CONGRESS
  ggnetwork::geom_edges(
    aes(size = connection_strength / max(connection_strength)), # Edge thickness proportional to connection strength, semi-transparent gray
    color="grey50", alpha = 0.3) +
  ggnetwork::geom_nodes(size = nodesize, 
                        aes(color = factor(color, levels = c("D", "R", "I")))) +
  ggplot2::theme(panel.grid = element_blank(),          # Remove grid lines
                 axis.text = element_blank(),           # Remove axis text
                 axis.ticks = element_blank(),          # Remove axis ticks
                 axis.title = element_blank(),          # Remove axis titles
                 panel.background = element_blank(),     # Remove panel background
                 strip.background = element_blank(),     # Remove facet backgrounds
                 strip.text = element_text(size = 10, face = "bold"),  # Style facet labels
                 legend.position = "bottom") +          # Position legend at bottom
  ggplot2::scale_size_continuous(guide = F, range = c(0.1, 2)) +
  ggplot2::scale_colour_manual(breaks = c("D", "R", "I"), 
                               values = c("#04668C", "#B9121B", "#588F27"), # Use standard political colors: Blue for Democrats, Red for Republicans, Green for Independents
                               labels = c("Democratic", "Republican", "Independent"),
                               name = "Party") +
  ggnetwork::geom_nodetext_repel(aes(label=text), size = 2.5, vjust = 0.7, hjust = -0.5, color = 'black')

# ENHANCE PLOT WITH SPECIFIC SENATOR LABELS
# Add custom labels for notable senators in specific Congressional sessions
pic <- basic_pic +
  geom_nodetext_repel(
    aes(label = vertex.names), color = 'black', 
    data = function(x) {
      x[x$vertex.names == "GILLIBRAND, Kirsten" &
          x$color == "D" &
          x$congress == "116-th Congress (2019-2021)", ]
    },
    vjust = -0.7, hjust = 1.0
  ) + 
  geom_nodetext_repel(
    aes(label = vertex.names), color = 'black',
    vjust = -2.6, hjust = 60, 
    size = 4, 
    data = function(x) {
      x[x$vertex.names == "HARRIS, Kamala Devi" &
          x$color == "D" &
          x$congress == "116-th Congress (2019-2021)",]
    }
  )

# DISPLAY FINAL VISUALIZATION
# Show the completed network plot with all enhancements
pic