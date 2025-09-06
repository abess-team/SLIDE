rm(list = ls()); gc(reset = TRUE)
path <- "~/voting/"
setwd(path)
library(ggnetwork)
library(ggpmisc)
library(network)
library(sna)
library(ggpubr)

# numbers: votes and senators
load("senate_model_data.rda")
tab_info <- sapply(senate_model_data, function(x) {
  dim(t(x[["vote"]]))
})
xtable::xtable(tab_info)

plot_data <- function(theta, dat, complete = FALSE) {
  act_set <- which(theta != 0, arr.ind = TRUE)
  act_set <- act_set[act_set[, 1] < act_set[, 2],]
  tmp <- act_set
  for(i in 1:nrow(act_set)) {
    for(j in 1:ncol(act_set)) {
      act_set[i, j] <- colnames(theta)[tmp[i, j]]
    }
  }
  act_set <- apply(act_set, 2, as.character)
  p <- ncol(theta)
  p1 <- length(unique(as.vector(act_set)))
  
  congress <- dat[["congress"]][1]
  
  nodesize <- 4
  
  if(complete) {
    g <- network::network.initialize(p, directed = FALSE)
    network::network.edgelist(tmp, g)
    network::network.vertex.names(g) <- colnames(theta)
  } else {
    g <- network::network(act_set, matrix.type = "edgelist", directed = FALSE)
  }
  
  connection_strength <- vector(length = nrow(act_set))
  for(i in 1:length(connection_strength)) {
    connection_strength[i] <- abs(theta[tmp[i, 1], tmp[i, 2]])
  }
  connection_strength <- as.numeric(connection_strength)
  g %v% "connection_strength" <- connection_strength
  
  vertex_name <- network::network.vertex.names(g)
  match_id <- match(vertex_name, dat$name)
  party <- dat$party[match_id]
  g %v% "color" <- as.character(party)
  
  if (congress == 116){
    g %v% "text" <- ifelse(g %v% "vertex.names" %in% c("HARRIS, Kamala Devi", "GILLIBRAND, Kirsten"), g %v% "vertex.names", "")
  } else {
    g %v% "text" <- ""
  }
  
  ggnetwork_default <- ggnetwork::ggnetwork(g, layout = "eigen")
  return(ggnetwork_default)
}
load("senate_nodewise.rda")
load("senate_model_data.rda")
index <- which(names(senate_nodewise) %in% as.character(congress_range))
senate_nodewise <- senate_nodewise[index]
senate_model_data <- senate_model_data[index]

pdat_list <- mapply(function(x, y) {
  y <- y[["info"]]
  colnames(x) <- y[["name"]]
  x
  pdat <- plot_data(x, y, complete = FALSE)
  start_year <- 2 * y[["congress"]][1] + 1787
  end_year <- start_year + 2
  congress_info <- paste0(sprintf("%s-th Congress", y[["congress"]][1]), 
                          " (", start_year, "-", end_year, ")")
  pdat <- cbind(pdat, "congress" = congress_info)
  pdat
}, senate_nodewise, senate_model_data, SIMPLIFY = FALSE)
pdat <- do.call("rbind.data.frame", pdat_list)

nodesize <- 3
set.seed(123)
basic_pic <- ggplot2::ggplot(pdat, aes(x, y, xend = xend, yend = yend)) +
  facet_wrap(congress ~ ., scales = "free", ncol = 2) + 
  ggnetwork::geom_edges(aes(size = connection_strength / max(connection_strength)), 
                        color="grey50", alpha = 0.3) +
  ggnetwork::geom_nodes(size = nodesize, 
                        aes(color = factor(color, levels = c("D", "R", "I")))) +
  ggplot2::theme(panel.grid = element_blank(),axis.text = element_blank(),
                 axis.ticks = element_blank(),axis.title = element_blank(),
                 panel.background = element_blank(),
                 strip.background = element_blank(), 
                 strip.text = element_text(size = 10, face = "bold"),
                 legend.position = "bottom") + 
  ggplot2::scale_size_continuous(guide = F, range = c(0.1, 2)) +
  ggplot2::scale_colour_manual(breaks = c("D", "R", "I"), 
                               values = c("#04668C", "#B9121B", "#588F27"),
                               labels = c("Democratic", "Republican", "Independent"),
                               name = "Party") +
  ggnetwork::geom_nodetext_repel(aes(label=text), size = 2.5, vjust = 0.7, hjust = -0.5, color = 'black')
pic <- basic_pic

republican_in_democratic <- subset(pdat, congress == "115-th Congress (2017-2019)", select = c(x, y, color, vertex.names))
republican_in_democratic <- subset(republican_in_democratic, x > median(x))
republican_in_democratic <- subset(republican_in_democratic, y > median(y))
republican_in_democratic <- unique(republican_in_democratic)
republican_in_democratic <- subset(republican_in_democratic, color == "R")
republican_in_democratic

democratic_in_republican <- subset(pdat, congress == "116-th Congress (2019-2021)", select = c(x, y, color, vertex.names))
democratic_in_republican <- subset(democratic_in_republican, x > median(x))
democratic_in_republican <- subset(democratic_in_republican, y > median(y))
democratic_in_republican <- unique(democratic_in_republican)
democratic_in_republican <- subset(democratic_in_republican, color == "D")
democratic_in_republican

pic <- basic_pic +
  geom_nodetext_repel(
    aes(label = vertex.names), color = 'black', 
    data = function(x) {
      x[x$vertex.names %in% c("FLAKE, Jeff", "HARRIS, Kamala Devi"), ]
    }, vjust = -0.7, hjust = 1.0
  )
pic

pic <- basic_pic +
  geom_nodetext_repel(
    aes(label = vertex.names), color = 'black', 
    data = function(x) {
      x[x$vertex.names == "FLAKE, Jeff" &
          x$color == "R" &
          x$congress == "115-th Congress (2017-2019)", ]
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
pic
