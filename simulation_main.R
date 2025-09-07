sim <- function(seed, type, n, p, method, alpha = 0.4, beta = NULL, degree = 3, thres = NULL) {
  runtime <- 0.0
  print(paste0("seed: ", seed))
  if (p <= 20) {
    train <- generate.bmn.data(n, p, type = type, seed = seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree)
    valid <- generate.bmn.data(n, theta = train[["theta"]], seed = 10000 + seed)
  } else {
    train <- generate.bmn.data(n, p, type = type, seed = seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree, method = "gibbs")
    valid <- generate.bmn.data(n, theta = train[["theta"]], seed = 10000 + seed, method = "gibbs")
  }
  sample_weight <- c(train[["weight"]], valid[["weight"]])
  theta <- train[["theta"]]
  edge_num_vec <- rowSums(theta != 0)
  train <- train[["data"]]
  valid <- valid[["data"]]
  
  if (save) {
    write.csv(theta, paste0(path, result_file, method, "_type", type, "_p", p, "_alpha", alpha, "_beta", beta, "_n", n, "_", mark, "_theta_true.csv"))
  }
  true_size <- 0.5 * sum(theta != 0)
  print(paste0("true size: ", true_size))
  
  cc <- 12
  shrink <- 0.8
  c_vector <- c()
  while (cc >= 1) {
    c_vector <- c(c_vector, cc)
    cc <- floor(cc * shrink)
  }
  
  if (is.null(thres)) {
    thres <- alpha / 2
  }
  c_RPLE <- 0.2
  c_RISE <- 0.4
  c_logRISE <- 0.8
  c_list <- 10^(seq(-4, 1, length.out = min(c(p - 2, 100))))
  train <- cbind(sample_weight[1:nrow(train)], train)
  valid <- cbind(sample_weight[(nrow(train) + 1):(nrow(train) + nrow(valid))], valid)
  
  if ("RPLE" %in% method)
    runtime <- system.time(
      RPLE <- ISO(train, valid, c_list = c_RPLE, method = "RPLE")
    )[3]
  else
    RPLE <- NULL
  if ("RPLE_thres" %in% method)
    runtime <- system.time(
      RPLE_thres <- ISO(train, valid, c_list = c_RPLE, thres = thres, method = "RPLE")
    )[3]
  else
    RPLE_thres <- NULL
  if ("RPLE_cv_thres" %in% method)
    runtime <- system.time(
      RPLE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "RPLE")
    )[3]
  else
    RPLE_cv_thres <- NULL
  
  if ("ELASSO" %in% method)
    runtime <- system.time(ELASSO <- ELASSO(train, valid))[3]
  else
    ELASSO <- NULL
  if ("ELASSO_thres" %in% method)
    runtime <- system.time(
      ELASSO_thres <- ELASSO(train, valid, thres = thres)
    )[3]
  else
    ELASSO_thres <- NULL
  
  if ("RISE" %in% method)
    runtime <- system.time(
      RISE <- ISO(train, valid, c_list = c_RISE, method = "RISE")
    )[3]
  else
    RISE <- NULL
  if ("RISE_thres" %in% method)
    runtime <- system.time(
      RISE_thres <- ISO(train, valid, c_list = c_RISE, thres = thres, method = "RISE")
    )[3]
  else
    RISE_thres <- NULL
  if ("RISE_cv_thres" %in% method)
    runtime <- system.time(
      RISE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "RISE")
    )[3]
  else
    RISE_cv_thres <- NULL
  
  if ("logRISE" %in% method) {
    runtime <- system.time(
      logRISE <- ISO(train, valid, c_list = c_logRISE, method = "logRISE")
    )[3]
  } else {
    logRISE <- NULL
  }
  if ("logRISE_thres" %in% method) {
    runtime <- system.time(
      logRISE_thres <- ISO(train, valid, c_list = c_logRISE, thres = thres, method = "logRISE")
    )[3]
  } else {
    logRISE_thres <- NULL
  }
  if ("logRISE_cv_thres" %in% method) {
    runtime <- system.time(
      logRISE_cv_thres <- ISO(train, valid, c_list = c_list, thres = thres, method = "logRISE")
    )[3]
  } else {
    logRISE_cv_thres <- NULL
  }
  
  if ("LogRelax" %in% method)
    runtime <- system.time(LogRelax <- LogRelax(train, valid))[3]
  else
    LogRelax <- NULL
  if ("LogRelaxTAlpha" %in% method)
    runtime <- system.time(LogRelaxTAlpha <- LogRelax(train, valid, theory_alpha =
                                                        TRUE))[3]
  else
    LogRelaxTAlpha <- NULL
  if ("LogRelax_thres" %in% method)
    runtime <- system.time(LogRelax_thres <- LogRelax(train, valid, thres = thres))[3]
  else
    LogRelax_thres <- NULL
  if ("LogRelaxTAlpha_thres" %in% method)
    runtime <- system.time(LogRelaxTAlpha_thres <- LogRelax(train, valid, theory_alpha = TRUE, thres = thres))[3]
  else
    LogRelaxTAlpha_thres <- NULL
  
  train <- train[, -1]
  valid <- valid[, -1]
  pool_data <- rbind(train, valid)
  if ("SLIDE_oracle" %in% method) {
    runtime <- system.time(
      SLIDE_oracle <- slide(pool_data, weight = sample_weight, max.support.size = edge_num_vec, tune.type = "gic", support.size = edge_num_vec)[[1]]
    )[3]
    diag(SLIDE_oracle) <- 0.0
  } else {
    SLIDE_oracle <- NULL
  }
  if ("SLIDE" %in% method) {
    runtime <- system.time(
      SLIDE <- slide(pool_data, weight = sample_weight, tune.type = "gic", ic.scale = 2, graph.threshold = thres)[[1]]
    )[3]
    diag(SLIDE) <- 0.0
  } else {
    SLIDE <- NULL
  }
  
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
  rate <- sim.rate(res, theta)
  print(rate)
  loss <- sim.loss(res, theta)
  is.rec <- sim.prop(res, theta)
  mse <- sim.mse(res, theta)
  if (thres == 0.0) {
    auc <- sim.auc(res, theta)
    res <- c(rate, loss, is.rec, mse, auc, runtime)
  } else {
    res <- c(rate, loss, is.rec, mse)
  }
  
  return(res)
}
