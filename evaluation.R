compute_rate <- function(X.opt, sparse) {
  X.opt <- as.matrix(X.opt)
  est_act <- which(X.opt != 0, arr.ind = TRUE)
  if(any(diag(X.opt) != 0)) {
    est_act <- est_act[-which(est_act[,1] == est_act[,2]),]
  }
  est_act <- as.data.frame(est_act)
  true_act <- which(sparse != 0, arr.ind = TRUE)
  if(any(diag(sparse) != 0)) {
    true_act <- true_act[-which(true_act[,1] == true_act[,2]),]
  }
  true_act <- as.data.frame(true_act)
  est_inact <- which(X.opt == 0, arr.ind = TRUE)
  if(any(diag(X.opt) == 0)) {
    est_inact <- est_inact[-which(est_inact[,1] == est_inact[,2]),]
  }
  est_inact <- as.data.frame(est_inact)
  true_inact <- which(sparse == 0, arr.ind = TRUE)
  if(any(diag(sparse) == 0)) {
    true_inact <- true_inact[-which(true_inact[,1] == true_inact[,2]),]
  }
  true_inact <- as.data.frame(true_inact)
  TP <- as.numeric(nrow(dplyr::intersect(est_act, true_act)))
  FP <- as.numeric(nrow(dplyr::intersect(est_act, true_inact)))
  TN <- as.numeric(nrow(dplyr::intersect(est_inact, true_inact)))
  FN <- as.numeric(nrow(dplyr::intersect(est_inact, true_act)))
  TPR <- TP / nrow(true_act)
  FPR <- FP / nrow(true_inact)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  return(c(TPR = TPR, FPR = FPR, MCC = MCC))
}

compute_loss <- function(X.opt, sigma_inv) {
  X.opt <- as.matrix(X.opt)
  loss_op <- norm(X.opt - sigma_inv, "2")
  loss_l1 <- norm(X.opt - sigma_inv, "1")
  loss_F <- norm(X.opt - sigma_inv, "F")
  return(c(loss_op = loss_op, loss_l1 = loss_l1, loss_F = loss_F))
}

compute_mse <- function(X.opt, sigma_inv) {
  X.opt <- as.matrix(X.opt)
  sigma_inv <- as.matrix(sigma_inv)
  temp <- X.opt - sigma_inv
  mse <- mean(temp[upper.tri(temp, diag = TRUE)]^2)
  return(mse = mse)
}

compute_auc <- function(X.opt, sigma_inv) {
  X.opt <- abs(as.matrix(X.opt))
  sigma_inv <- as.matrix(sigma_inv)
  sigma_inv[abs(sigma_inv) > 0.0] <- 1
  X.opt <- X.opt[upper.tri(X.opt)]
  sigma_inv <- sigma_inv[upper.tri(sigma_inv)]
  auc = pROC::roc(sigma_inv, X.opt)[["auc"]][1]
  return(auc = auc)
}

sim.rate <- function(res, sparse){
  sparse <- as.matrix(sparse)
  rate <- c()
  for(i in 1:length(res)){
    if(!is.null(res[[i]])){
      temp <- as.matrix(res[[i]])
      temp_rate <- compute_rate(temp, sparse)
      names(temp_rate) <- paste(names(res)[i], names(temp_rate), sep = "_")
      rate <- c(rate, temp_rate)
    }
  }
  return(rate)
}

sim.prop <- function(res, sparse){
  sparse <- as.matrix(sparse)
  is.recov <- c()
  for(i in 1: length(res)){
    if(!is.null(res[[i]])){
      temp <- as.matrix(res[[i]])
      temp_rec <- identical(which(temp != 0), which(sparse != 0))
      names(temp_rec) <- paste(names(res)[i], names(temp_rec), sep = "_")
      is.recov <- c(is.recov, temp_rec)
    }
  }
  return(is.recov)
}

sim.loss <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  loss <- c()
  for(i in 1:length(res)) {
    if(!is.null(res[[i]])) {
      temp <- as.matrix(res[[i]])
      temp_loss <- compute_loss(temp, sigma_inv)
      names(temp_loss) <- paste(names(res)[i], names(temp_loss), sep = "_")
      loss <- c(loss, temp_loss)
    }
  }
  return(loss)
}

sim.mse <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  mse <- c()
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

sim.auc <- function(res, sigma_inv) {
  sigma_inv <- as.matrix(sigma_inv)
  auc <- c()
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