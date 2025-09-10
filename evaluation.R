# ======================================================================================
# EVALUATION FUNCTIONS FOR ISING MODEL RECONSTRUCTION ASSESSMENT
# ======================================================================================
# This file contains functions to evaluate the performance of Ising model reconstruction
# methods by computing various metrics like true positive rate, false positive rate, and
# F norm.


#' Compute Graphical Structure Recovery Performance
#' 
#' Calculates True Positive Rate (TPR), False Positive Rate (FPR), and 
#' Matthews Correlation Coefficient (MCC) for binary edge classification
#' in graphical models.
#' 
#' @param X.opt Estimated precision matrix (or adjacency matrix)
#' @param sparse True sparse precision matrix (ground truth)
#' @return Named vector with TPR, FPR, and MCC values
#' 
#' TPR: True Positive Rate = TP / (TP + FN)
#' FPR: False Positive Rate = FP / (FP + TN) 
#' MCC: Matthews Correlation Coefficient = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
compute_rate <- function(X.opt, sparse) {
  # Convert inputs to matrices for consistent handling
  X.opt <- as.matrix(X.opt)

  # Find estimated active edges (non-zero elements), excluding diagonal
  est_act <- which(X.opt != 0, arr.ind = TRUE)
  if(any(diag(X.opt) != 0)) {
    est_act <- est_act[-which(est_act[,1] == est_act[,2]),]
  }
  est_act <- as.data.frame(est_act)
  
  # Find true active edges (non-zero elements in ground truth), excluding diagonal
  true_act <- which(sparse != 0, arr.ind = TRUE)
  if(any(diag(sparse) != 0)) {
    true_act <- true_act[-which(true_act[,1] == true_act[,2]),]
  }
  true_act <- as.data.frame(true_act)

  # Find estimated inactive edges (zero elements), excluding diagonal
  est_inact <- which(X.opt == 0, arr.ind = TRUE)
  if(any(diag(X.opt) == 0)) {
    est_inact <- est_inact[-which(est_inact[,1] == est_inact[,2]),]
  }
  est_inact <- as.data.frame(est_inact)
  
  # Find true inactive edges (zero elements in ground truth), excluding diagonal
  true_inact <- which(sparse == 0, arr.ind = TRUE)
  if(any(diag(sparse) == 0)) {
    true_inact <- true_inact[-which(true_inact[,1] == true_inact[,2]),]
  }
  true_inact <- as.data.frame(true_inact)

  # Calculate confusion matrix components
  TP <- as.numeric(nrow(dplyr::intersect(est_act, true_act)))     # True Positives
  FP <- as.numeric(nrow(dplyr::intersect(est_act, true_inact)))   # False Positives
  TN <- as.numeric(nrow(dplyr::intersect(est_inact, true_inact))) # True Negatives
  FN <- as.numeric(nrow(dplyr::intersect(est_inact, true_act)))   # False Negatives
  
  # Calculate performance metrics
  TPR <- TP / nrow(true_act)
  FPR <- FP / nrow(true_inact)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  return(c(TPR = TPR, FPR = FPR, MCC = MCC))
}

#' Compute Loss Between Estimated and True Matrices
#' 
#' Calculates various matrix norms to measure the difference between
#' estimated and true precision matrices.
#' 
#' @param X.opt Estimated precision matrix
#' @param sigma_inv True precision matrix (inverse covariance)
#' @return Named vector with operator norm, L1 norm, and Frobenius norm losses
compute_loss <- function(X.opt, sigma_inv) {
  X.opt <- as.matrix(X.opt)

  # Calculate different matrix norms as loss measures
  loss_op <- norm(X.opt - sigma_inv, "2")  # Operator norm (largest singular value)
  loss_l1 <- norm(X.opt - sigma_inv, "1")  # L1 norm (sum of absolute values)
  loss_F <- norm(X.opt - sigma_inv, "F")   # Frobenius norm (sqrt of sum of squares)

  return(c(loss_op = loss_op, loss_l1 = loss_l1, loss_F = loss_F))
}


#' Compute Mean Squared Error for Upper Triangular Elements
#' 
#' Calculates MSE considering only the upper triangular part of the matrices
#' (including diagonal), since precision matrices are symmetric.
#' 
#' @param X.opt Estimated precision matrix
#' @param sigma_inv True precision matrix
#' @return Scalar MSE value
compute_mse <- function(X.opt, sigma_inv) {
  X.opt <- as.matrix(X.opt)
  sigma_inv <- as.matrix(sigma_inv)

  # Calculate element-wise difference
  temp <- X.opt - sigma_inv

  # Compute MSE using only upper triangular elements (including diagonal)
  # This avoids double-counting due to symmetry
  mse <- mean(temp[upper.tri(temp, diag = TRUE)]^2)
  return(mse = mse)
}

#' Compute Area Under the ROC Curve for Edge Detection
#' 
#' Treats edge detection as a binary classification problem and computes AUC
#' using the absolute values of estimated coefficients as prediction scores.
#' 
#' @param X.opt Estimated precision matrix (used as prediction scores)
#' @param sigma_inv True precision matrix (used as binary labels)
#' @return AUC value between 0 and 1
compute_auc <- function(X.opt, sigma_inv) {
  X.opt <- abs(as.matrix(X.opt))
  sigma_inv <- as.matrix(sigma_inv)

  # Convert true matrix to binary labels (0 = no edge, 1 = edge)
  sigma_inv[abs(sigma_inv) > 0.0] <- 1

  # Extract upper triangular elements to avoid symmetry duplication
  X.opt <- X.opt[upper.tri(X.opt)]
  sigma_inv <- sigma_inv[upper.tri(sigma_inv)]

  # Compute ROC curve and extract AUC
  auc = pROC::roc(sigma_inv, X.opt)[["auc"]][1]
  
  return(auc = auc)
}

# ==============================================================================
# WRAPPER FUNCTIONS FOR MULTIPLE METHOD COMPARISON
# ==============================================================================
# These functions below apply the above evaluation metrics to a list of results
# from different estimation methods, enabling easy comparison.

#' Compute TPR, FPR and MCC to Multiple Methods
#' 
#' Wrapper function that computes TPR, FPR, and MCC for each method
#' in a list of estimation results.
#' 
#' @param res List of estimated precision matrices from different methods
#' @param sparse True sparse precision matrix (ground truth)
#' @return Named vector with rates for each method (method_name_metric_name)
sim.rate <- function(res, sparse){
  sparse <- as.matrix(sparse)
  rate <- c()

  # Iterate through each method's results
  for(i in 1:length(res)){
    if(!is.null(res[[i]])){
      temp <- as.matrix(res[[i]])
      # Compute rates for this method
      temp_rate <- compute_rate(temp, sparse)
      # Create descriptive names combining method and metric names
      names(temp_rate) <- paste(names(res)[i], names(temp_rate), sep = "_")
      rate <- c(rate, temp_rate)
    }
  }
  return(rate)
}

#' Check Exact Graphical Structure Recovery for Multiple Methods
#' 
#' Determines whether each method exactly recovers the true sparsity pattern
#' (i.e., identifies exactly the same set of non-zero elements).
#' 
#' @param res List of estimated precision matrices from different methods
#' @param sparse True sparse precision matrix (ground truth)
#' @return Named logical vector indicating exact recovery for each method
sim.prop <- function(res, sparse){
  sparse <- as.matrix(sparse)
  is.recov <- c()

  # Check exact sparsity pattern recovery for each method
  for(i in 1: length(res)){
    if(!is.null(res[[i]])){
      temp <- as.matrix(res[[i]])
      # Check if non-zero positions exactly match the true pattern
      temp_rec <- identical(which(temp != 0), which(sparse != 0))
      names(temp_rec) <- paste(names(res)[i], names(temp_rec), sep = "_")
      is.recov <- c(is.recov, temp_rec)
    }
  }
  return(is.recov)
}

#' Apply Loss Function Computation to Multiple Methods
#' 
#' Wrapper function that computes operator, L1, and Frobenius norm losses
#' for each method in a list of estimation results.
#' 
#' @param res List of estimated precision matrices from different methods
#' @param sigma_inv True precision matrix (ground truth)
#' @return Named vector with losses for each method (method_name_loss_type)
sim.loss <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  loss <- c()

  # Compute losses for each method
  for(i in 1:length(res)) {
    if(!is.null(res[[i]])) {
      temp <- as.matrix(res[[i]])
      temp_loss <- compute_loss(temp, sigma_inv)
      # Create descriptive names combining method and loss type
      names(temp_loss) <- paste(names(res)[i], names(temp_loss), sep = "_")
      loss <- c(loss, temp_loss)
    }
  }
  return(loss)
}

#' Apply MSE Computation to Multiple Methods
#' 
#' Wrapper function that computes mean squared error for each method
#' in a list of estimation results.
#' 
#' @param res List of estimated precision matrices from different methods
#' @param sigma_inv True precision matrix (ground truth)
#' @return Named vector with MSE for each method
sim.mse <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  mse <- c()

  # Compute MSE for each method
  for(i in 1:length(res)) {
    if(!is.null(res[[i]])) {
      temp <- as.matrix(res[[i]])
      temp_mse <- compute_mse(temp, sigma_inv)
      names(temp_mse) <- paste(names(res[i]), "mse", sep = "_")
      mse <- c(mse, temp_mse)
    }
  }
  return(mse)
}

#' Apply AUC Computation to Multiple Methods
#' 
#' Wrapper function that computes area under the ROC curve for each method
#' in a list of estimation results.
#' 
#' @param res List of estimated precision matrices from different methods
#' @param sigma_inv True precision matrix (ground truth)
#' @return Named vector with AUC for each method
sim.auc <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  auc <- c()
  
  # Compute AUC for each method
  for(i in 1:length(res)) {
    if(!is.null(res[[i]])) {
      temp <- as.matrix(res[[i]])
      temp_auc <- compute_auc(temp, sigma_inv)
      names(temp_auc) <- paste(names(res[i]), "auc", sep = "_")
      auc <- c(auc, temp_auc)
    }
  }
  return(auc)
}