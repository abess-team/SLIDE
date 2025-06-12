# Optimal Size
rm(list = ls()); gc(reset = TRUE)
library(reshape2)
library(ggplot2)
library(ggpmisc)

# path <- "D:\\ising-L0"
# path <- "~/splicing-ising/code-simulate/0516"
path <- "/Users/zhujin/splicing-ising/code-simulate/0910"
# path <- "C:/Users/ZHUJ68/SLIDE/"
setwd(path)

# method_list <- c("RPLE_thres", "nodewise_logistic_gic2")
# method_list <- c("RPLE_thres")
# method_list <- c("RPLE_thres", "logRISE_thres")
method_list <- c("RPLE_thres", "RISE_thres", "logRISE_thres", "nodewise_logistic_gic2")

# case_list <- 2
case_list <- c(1, 2)
res <- list()
for (case in case_list) {
  if (case == 1) {
    type <- 10

    degree_list <- 4
    
    alpha_list <- 0.3
    beta_list <- 0.4
    p_list <- c(9, 16, 25, 36, 49, 64, 81, 100)
  } else {
    type <- 8
    
    degree_list <- c(3)
    
    alpha_list <- 0.3
    beta_list <- 0.6
    
    p_list <- c(8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52)
  }
  
  n_matrix <- matrix(0, length(method_list), length(p_list))
  for(i in 1: length(method_list)) {
    for(j in 1:length(p_list)) {
      name <- paste0("result_p/", method_list[i], "_", "type", type, 
                     "_p", p_list[j], "_alpha", alpha_list, "_beta", beta_list, 
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
  
  res[[as.character(case)]] <- data.frame(t(n_matrix), "p" = p_list, "type" = type)
}
pdat <- do.call("rbind.data.frame", res)

method_name <- method_list
method_name <- sapply(strsplit(method_name, split = "_"), function(x) x[[1]])
method_name[method_name == "nodewise"] <- "SLIDE"


if (length(case_list) > 1) {
  colnames(pdat)[1:length(method_name)] <- method_name
  data_melt <- melt(pdat, id.vars = c("type", "p"), 
                    variable.name = "Method", value.name = "n")
  data_melt[["type"]] <- ifelse(data_melt[["type"]] == 8, 
                                "Random regular graph", 
                                "Square lattice")
  
  formula <- y ~ log(x)
  pj <- position_jitter(width = 1, height = 0)
  p <- ggplot(data = subset(data_melt), aes(x = p, y = n, group = Method)) + 
    facet_wrap(type ~ ., ncol = 2, scales = "free") + 
    geom_point(aes(color = Method, shape = Method), size = 2, position = pj) +
    # geom_smooth(aes(color = Method), se = FALSE, 
    #             formula = formula, method = "lm") +
    stat_poly_line(aes(color = Method), formula = formula, se = FALSE) +
    stat_poly_eq(aes(color = Method, label = after_stat(eq.label)),
                 formula = formula) +
    xlab(expression(italic(p)))+
    ylab("Optimal sample size")+
    scale_color_discrete(name = "Methods") +
    scale_shape_discrete(name = "Methods") + 
    theme_bw() + 
    theme(legend.position = "bottom", 
          legend.box.margin = margin(-12, 0, -10, 0))
  p
  ggsave(p, filename = sprintf("p_sc.png", type), width = 8, height = 3)  
  ggsave(p, filename = sprintf("p_sc.pdf", type), width = 8, height = 3)  
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
  
  ggsave(p, filename = sprintf("p_type%s.pdf", type), width = 8, height = 4)  
}


