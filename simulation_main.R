sim <- function(seed,
                type,
                n,
                p,
                method,
                alpha = 0.4,
                beta = NULL,
                degree = 3,
                thres = NULL) {
  runtime <- 0.0
  print(paste0("seed: ", seed))
  if (p <= 20) {
    # train <- generate.bmn.freq.data(n, p, type = type, seed = seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree)
    # valid <- generate.bmn.freq.data(n, p, type = type, seed = 10000 * seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree)
    train <- generate.bmn.data(
      n,
      p,
      type = type,
      seed = seed,
      graph.seed = 1,
      alpha = alpha,
      beta = beta,
      degree = degree
    )
    # valid <- generate.bmn.data(n, p, type = type, seed = 10000 * seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree)
    valid <- generate.bmn.data(n, theta = train[["theta"]], seed = 10000 + seed)
  } else {
    train <- generate.bmn.data(
      n,
      p,
      type = type,
      seed = seed,
      graph.seed = 1,
      alpha = alpha,
      beta = beta,
      degree = degree,
      method = "gibbs"
    )
    # valid <- generate.bmn.data(n, p, type = type, seed = 10000 * seed, graph.seed = 1, alpha = alpha, beta = beta, degree = degree, method = "gibbs")
    valid <- generate.bmn.data(n,
                               theta = train[["theta"]],
                               seed = 10000 + seed,
                               method = "gibbs")
  }
  sample_weight <- c(train[["weight"]], valid[["weight"]])
  theta <- train[["theta"]]
  edge_num_vec <- rowSums(theta != 0)
  train <- train[["data"]]
  valid <- valid[["data"]]
  
  if (save) {
    write.csv(
      theta,
      paste0(
        path,
        result_file,
        method,
        "_type",
        type,
        "_p",
        p,
        "_alpha",
        alpha,
        "_beta",
        beta,
        "_n",
        n,
        "_",
        mark,
        "_theta_true.csv"
      )
    )
  }
  true_size <- 0.5 * sum(theta != 0)
  print(paste0("true size: ", true_size))
  # c_vector <- c(12, 10, 8, 6, 4, 2, 1)
  # c_vector <- c(8, 4, 2, 1)
  # c_vector <- c(12, 9, 7, 5, 4, 3, 2, 1)
  
  cc <- 12
  shrink <- 0.8
  
  #cc <- 5; shrink <- 0.5;
  #cc <- round(true_size / 2) - 1; shrink <- 0.8;
  #cc <- 16; shrink <- 0.5;
  c_vector <- c()
  while (cc >= 1) {
    c_vector <- c(c_vector, cc)
    cc <- floor(cc * shrink)
  }
  
  # c_type = 1的效果更好
  if (is.null(thres)) {
    thres <- alpha / 2
  }
  cluster <- 3
  if ("splicing" %in% method) {
    runtime <- system.time(
      splicing <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 1,
        cv = FALSE,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing <- NULL
  }
  if ("splicing2" %in% method) {
    runtime <- system.time(
      splicing2 <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 2,
        cv = FALSE,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing2 <- NULL
  }
  if ("splicing_kmeans" %in% method) {
    runtime <- system.time(
      splicing_kmeans <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 1,
        cv = FALSE,
        cluster = cluster,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing_kmeans <- NULL
  }
  if ("splicing_thres" %in% method) {
    runtime <- system.time(
      splicing_thres <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 1,
        thres = thres,
        cv = FALSE,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing_thres <- NULL
  }
  
  if ("splicing_cv" %in% method) {
    runtime <- system.time(
      splicing_cv <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 1,
        cv = TRUE,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing_cv <- NULL
  }
  if ("splicing_cv_thres" %in% method) {
    runtime <- system.time(
      splicing_cv_thres <- splicing_whole(
        train,
        valid,
        type = type,
        true_size = true_size,
        c_vector = c_vector,
        c_type = 1,
        thres = thres,
        cv = TRUE,
        warmstart = TRUE,
        opt_method = 3,
        seed = seed
      )
    )[3]
  } else {
    splicing_cv_thres <- NULL
  }
  
  c_RPLE <- 0.2
  c_RISE <- 0.4
  c_logRISE <- 0.8
  # c_list <- c(0:5) * 0.2
  c_list <- 10 ^ (seq(-4, 1, length.out = min(c(p - 2, 100))))
  train <- cbind(sample_weight[1:nrow(train)], train)
  valid <- cbind(sample_weight[(nrow(train) + 1):(nrow(train) + nrow(valid))], valid)
  
  if ("RPLE" %in% method)
    runtime <- system.time(RPLE <- ISO(train, valid, c_list = c_RPLE, method = "RPLE"))[3]
  else
    RPLE <- NULL
  if ("RPLE_kmeans" %in% method)
    runtime <- system.time(RPLE_kmeans <- ISO(
      train,
      valid,
      c_list = c_RPLE,
      cluster = cluster,
      method = "RPLE"
    ))[3]
  else
    RPLE_kmeans <- NULL
  if ("RPLE_thres" %in% method)
    runtime <- system.time(RPLE_thres <- ISO(
      train,
      valid,
      c_list = c_RPLE,
      thres = thres,
      method = "RPLE"
    ))[3]
  else
    RPLE_thres <- NULL
  if ("RPLE_cv_thres" %in% method)
    runtime <- system.time(RPLE_cv_thres <- ISO(
      train,
      valid,
      c_list = c_list,
      thres = thres,
      method = "RPLE"
    ))[3]
  else
    RPLE_cv_thres <- NULL
  
  if ("ELASSO" %in% method)
    runtime <- system.time(ELASSO <- ELASSO(train, valid))[3]
  else
    ELASSO <- NULL
  if ("ELASSO_thres" %in% method)
    runtime <- system.time(ELASSO_thres <- ELASSO(train, valid, thres = thres))[3]
  else
    ELASSO_thres <- NULL
  
  if ("RISE" %in% method)
    runtime <- system.time(RISE <- ISO(train, valid, c_list = c_RISE, method = "RISE"))[3]
  else
    RISE <- NULL
  if ("RISE_kmeans" %in% method)
    runtime <- system.time(RISE_kmeans <- ISO(
      train,
      valid,
      c_list = c_RISE,
      cluster = cluster,
      method = "RISE"
    ))[3]
  else
    RISE_kmeans <- NULL
  if ("RISE_thres" %in% method)
    runtime <- system.time(RISE_thres <- ISO(
      train,
      valid,
      c_list = c_RISE,
      thres = thres,
      method = "RISE"
    ))[3]
  else
    RISE_thres <- NULL
  if ("RISE_cv_thres" %in% method)
    runtime <- system.time(RISE_cv_thres <- ISO(
      train,
      valid,
      c_list = c_list,
      thres = thres,
      method = "RISE"
    ))[3]
  else
    RISE_cv_thres <- NULL
  
  if ("logRISE" %in% method) {
    runtime <- system.time(logRISE <- ISO(train, valid, c_list = c_logRISE, method = "logRISE"))[3]
  } else {
    logRISE <- NULL
  }
  if ("logRISE_kmeans" %in% method) {
    runtime <- system.time(
      logRISE_kmeans <- ISO(
        train,
        valid,
        c_list = c_logRISE,
        cluster = cluster,
        method = "logRISE"
      )
    )[3]
  } else {
    logRISE_kmeans <- NULL
  }
  if ("logRISE_thres" %in% method) {
    runtime <- system.time(
      logRISE_thres <- ISO(
        train,
        valid,
        c_list = c_logRISE,
        thres = thres,
        method = "logRISE"
      )
    )[3]
  } else {
    logRISE_thres <- NULL
  }
  if ("logRISE_cv_thres" %in% method) {
    runtime <- system.time(
      logRISE_cv_thres <- ISO(
        train,
        valid,
        c_list = c_list,
        thres = thres,
        method = "logRISE"
      )
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
    runtime <- system.time(LogRelaxTAlpha_thres <- LogRelax(train, valid, theory_alpha =
                                                              TRUE, thres = thres))[3]
  else
    LogRelaxTAlpha_thres <- NULL
  
  ## new implementation:
  # train <- train[train[, 1] > 0, ]
  # valid <- valid[valid[, 1] > 0, ]
  train <- train[, -1]
  valid <- valid[, -1]
  foldid <- c(rep(1, nrow(train)), rep(2, nrow(valid)))
  pool_data <- rbind(train, valid)
  # sample_weight <- pool_data[, 1]
  if ("splicing3" %in% method) {
    runtime <- system.time(
      splicing3 <- abessbmn(
        pool_data,
        weight = sample_weight,
        tune.type = "gic",
        support.size = (true_size - 10):true_size,
        c.max = round(2 * p / 3)
      )[["omega"]]
    )[3]
    splicing3 <- splicing3[, , dim(splicing3)]
  } else {
    splicing3 <- NULL
  }
  if ("splicing_cv_thres2" %in% method) {
    runtime <- system.time(
      splicing_cv_thres2 <- abessbmn(
        pool_data,
        weight = sample_weight,
        tune.type = "cv",
        nfolds = 2,
        foldid = foldid,
        c.max = round(2 * p / 3),
        graph.threshold = thres
      )[["optimal.omega"]]
    )[3]
  } else {
    splicing_cv_thres2 <- NULL
  }
  # if ("nodewise_logistic" %in% method) splicing_logistic <- nodewise_L0(pool_data,
  #                                                                       max.support.size = edge_num_vec,
  #                                                                       tune.type = "gic", support.size = edge_num_vec)[[1]] else splicing_logistic <- NULL
  # if ("nodewise_logistic_cv" %in% method) splicing_cv_logistic <- nodewise_L0(pool_data, tune.type = "cv",
  #                                                                             foldid = foldid, graph.threshold = thres)[[1]] else splicing_cv_logistic <- NULL
  if ("nodewise_logistic_oraclesize" %in% method) {
    runtime <- system.time(
      splicing_logistic <- nodewise_L0(
        pool_data,
        weight = sample_weight,
        max.support.size = edge_num_vec,
        tune.type = "gic",
        support.size = edge_num_vec
      )[[1]]
    )[3]
    diag(splicing_logistic) <- 0.0
  } else {
    splicing_logistic <- NULL
  }
  if ("nodewise_logistic_gic" %in% method) {
    runtime <- system.time(
      splicing_gic_logistic <- nodewise_L0(
        pool_data,
        weight = sample_weight,
        tune.type = "gic",
        ic.scale = 2
      )[[1]]
    )[3]
    diag(splicing_gic_logistic) <- 0.0
  } else {
    splicing_gic_logistic <- NULL
  }
  if ("nodewise_logistic_gic2" %in% method) {
    runtime <- system.time(
      splicing_gic_logistic <- nodewise_L0(
        pool_data,
        weight = sample_weight,
        tune.type = "gic",
        ic.scale = 1,
        graph.threshold = thres
      )[[1]]
    )[3]
    diag(splicing_gic_logistic) <- 0.0
  } else {
    splicing_gic_logistic <- NULL
  }
  if ("nodewise_logistic_bic" %in% method) {
    runtime <- system.time(
      splicing_bic_logistic <- nodewise_L0(
        pool_data,
        weight = sample_weight,
        tune.type = "bic",
        ic.scale = 1,
        graph.threshold = thres
      )[[1]]
    )[3]
    diag(splicing_bic_logistic) <- 0.0
  } else {
    splicing_bic_logistic <- NULL
  }
  if ("nodewise_logistic_cv" %in% method) {
    runtime <- system.time(
      splicing_cv_logistic <- nodewise_L0(
        pool_data,
        weight = sample_weight,
        tune.type = "cv",
        foldid = foldid,
        graph.threshold = thres
      )[[1]]
    )[3]
    diag(splicing_cv_logistic) <- 0.0
  } else {
    splicing_cv_logistic <- NULL
  }
  
  res <- list(
    splicing = splicing,
    splicing2 = splicing2,
    splicing_kmeans = splicing_kmeans,
    splicing_thres = splicing_thres,
    splicing_cv = splicing_cv,
    splicing_cv_thres = splicing_cv_thres,
    RPLE = RPLE,
    RPLE_kmeans = RPLE_kmeans,
    RPLE_thres = RPLE_thres,
    RPLE_cv_thres = RPLE_cv_thres,
    ELASSO = ELASSO,
    ELASSO_thres = ELASSO_thres,
    RISE = RISE,
    RISE_kmeans = RISE_kmeans,
    RISE_thres = RISE_thres,
    RISE_cv_thres = RISE_cv_thres,
    logRISE = logRISE,
    logRISE_kmeans = logRISE_kmeans,
    logRISE_thres = logRISE_thres,
    logRISE_cv_thres = logRISE_cv_thres,
    LogRelax = LogRelax,
    LogRelaxTAlpha = LogRelaxTAlpha,
    LogRelax_thres = LogRelax_thres,
    LogRelaxTAlpha_thres = LogRelaxTAlpha_thres,
    splicing3 = splicing3,
    splicing_cv_thres2 = splicing_cv_thres2,
    splicing_logistic = splicing_logistic,
    splicing_gic_logistic = splicing_gic_logistic,
    splicing_bic_logistic = splicing_bic_logistic,
    splicing_cv_logistic = splicing_cv_logistic
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
