# Optimal Size
rm(list = ls()); gc(reset = TRUE)
library(reshape2)
library(ggplot2)
library(ggpmisc)

# path <- "D:\\ising-L0"
# path <- "~/splicing-ising/code-simulate/0516"
path <- "/Users/zhujin/splicing-ising/code-simulate/0910/result_ws"
# path <- "C:/Users/ZHUJ68/SLIDE/result_ws/"
setwd(path)

# method_list <- c("nodewise_logistic_gic2")
# method_list <- c("RPLE_thres")
# method_list <- c("RPLE_thres", "logRISE_thres")
method_list <- c("RPLE_thres", "ELASSO_thres",
                 "RISE_thres", "logRISE_thres", "nodewise_logistic_gic2")

p_list <- c(16)
# case_list <- 1
case_list <- c(1, 2)
res <- list()
for (case in case_list) {
  if (case == 1) {
    omega <- 2.0
    degree_list <- c(3)
    alpha_list <- seq(30, 9, length=8) / 50
    beta_list <- (omega - alpha_list) / (degree_list - 1)
    type <- 8
  } else {
    omega <- 0.9
    degree_list <- c(4)
    alpha_list <-3 * (seq(30, 3, length=10) / 400)
    beta_list <- (omega - alpha_list) / (degree_list - 1)
    type <- 10
  }
  
  n_matrix <- matrix(0, length(method_list), length(alpha_list))
  for(i in 1: length(method_list)) {
    for(j in 1:length(alpha_list)) {
      name <- paste0(method_list[i], "_", 
                     "type", type, "_p", p_list, 
                     "_alpha", alpha_list[j], "_beta", beta_list[j], 
                     "_degree", degree_list, "_test.csv")
      
      tmp <- read.csv(name)
      temp <- colnames(tmp)[- 1]
      
      n_list <- sapply(strsplit(temp, "_"), 
                       function(x) {x[2]})
      n_list <- as.numeric(n_list)
      # n_max <- max(n_list)
      n_max <- n_list[length(n_list)]
      n_matrix[i, j] <- n_max
    }
  }
  
  res[[as.character(case)]] <- data.frame(t(n_matrix), 
                                          "alpha" = alpha_list, "type" = type)
}
pdat <- do.call("rbind.data.frame", res)

method_name <- method_list
method_name <- sapply(strsplit(method_name, split = "_"), function(x) x[[1]])
method_name[method_name == "nodewise"] <- "SLIDE"


if (length(case_list) > 1) {
  colnames(pdat)[1:length(method_name)] <- method_name
  data_melt <- melt(pdat, id.vars = c("type", "alpha"), 
                    variable.name = "Method", value.name = "n")
  data_melt[["type"]] <- ifelse(data_melt[["type"]] == 8, 
                                "Random regular graph", 
                                "Square lattice")
  formula <- y ~ x + I(x^2)
  p <- ggplot(data = data_melt, aes(x = 1 / alpha, y = n, group = Method)) + 
    facet_wrap(type ~ ., ncol = 2, scales = "free") + 
    geom_point(aes(color = Method, shape = Method), size = 2)+
    stat_poly_line(aes(color = Method), formula = formula, se = FALSE) + 
    stat_poly_eq(aes(color = Method, label = after_stat(eq.label)), 
                 formula = formula) + 
    xlab(expression(italic(lambda^{-1})))+
    ylab("Optimal sample size")+
    scale_color_discrete(name = "Methods") +
    scale_shape_discrete(name = "Methods") + 
    theme_bw() + 
    theme(legend.position = "bottom", 
          legend.box.margin = margin(-12, 0, -10, 0))
  p
  ggsave(p, filename = sprintf("ws.png", type), width = 8, height = 4)  
  ggsave(p, filename = sprintf("ws.pdf", type), width = 8, height = 4)  
} else {
  data <- data.frame(t(n_matrix))
  colnames(data) <- method_name
  
  data[["p"]] <- p_list
  data_melt <- melt(data, id.vars = "p", variable.name = "Method", value.name = "n")
  
  p <- ggplot(data = data_melt, aes(x = p, y = n, group = Method)) + 
    # geom_line(aes(color = Method), size = 0.5)+
    geom_point(aes(color = Method, shape = Method), size = 2)+
    geom_smooth(aes(color = Method), se = FALSE, formula = y ~ log(x), method = "lm") + 
    xlab("degree")+
    ylab("Optimal sample size")+
    scale_color_discrete(name = "methods")+
    scale_shape_discrete(name = "methods")
  p
  
  ggsave(p, filename = sprintf("ws_type%s.pdf", type), width = 8, height = 4)  
}


