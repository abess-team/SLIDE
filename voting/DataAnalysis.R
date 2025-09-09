rm(list = ls()); gc(reset = TRUE)
path <- here::here(); setwd(path)

library(reshape2)
library(abess)
library(dplyr)
library(ggnetwork)
library(ggpmisc)
library(network)
library(sna)
library(ggpubr)

congress_study <- 112:117
year <- seq(2011, 2023, by = 2)

senate_info <- read.csv("Sall_members.csv")
senate_info <- senate_info[senate_info[["congress"]] %in% congress_study, ]
senate_info[["party_code"]][!(senate_info[["party_code"]] %in% c(100, 200))] <- 300
senate_info[["party"]][senate_info[["party_code"]] == 100] <- "D"
senate_info[["party"]][senate_info[["party_code"]] == 200] <- "R"
senate_info[["party"]][senate_info[["party_code"]] == 300] <- "I"
senate_info[["year"]] <- 2 * senate_info[["congress"]] + 1789
senate_info[["age"]] <- senate_info[["year"]] - senate_info[["born"]]
table(senate_info[["party"]])
senate_info <- senate_info[c("congress", "chamber", "icpsr", "bioname", "party", "age")]
senate_info <- split(senate_info, f = senate_info[["congress"]])

votes <- read.csv("Sall_votes.csv")
votes <- votes[votes[["congress"]] %in% congress_study, ]
unique(votes[["chamber"]])
votes[["chamber"]] <- NULL
votes[["prob"]] <- NULL
table(votes[["cast_code"]])
votes[["cast_code"]] <- ifelse(votes[["cast_code"]] == 1, 1, -1)
votes <- split(votes, f = votes[["congress"]])
votes <- lapply(votes, function(x) {
  vote_x <- dcast(x, formula = icpsr ~ rollnumber, value.var = "cast_code")
  rownames(vote_x) <- vote_x[["icpsr"]]
  vote_x[["icpsr"]] <- NULL
  
  exclude_index <- which(rowMeans(is.na(vote_x)) > 0.0)
  if (length(exclude_index) != 0) {
    vote_x <- vote_x[-exclude_index, ]
  }
  
  exclude_index <- which(apply(vote_x, 2, function(x) {
    length(unique(x))
  }) == 1)
  if (length(exclude_index) != 0) {
    vote_x <- vote_x[, -exclude_index]
  }
  vote_x
})
sapply(votes, anyNA)
sapply(votes, dim)

senate_model_data <- mapply(function(x, y) {
  index <- match(as.numeric(rownames(y)), x[["icpsr"]])
  x <- x[index, ]
  x[["icpsr"]] <- NULL
  colnames(x)[colnames(x) == "bioname"] <- "name"
  list("info" = x, "vote" = as.matrix(y))
}, senate_info, votes, SIMPLIFY = FALSE)
save(senate_model_data, file = "senate_model_data.rda")

# numbers: votes and senators
tab_info <- sapply(senate_model_data, function(x) {
  dim(t(x[["vote"]]))
})
print(tab_info)

load("senate_model_data.rda")
senate_nodewise <- lapply(senate_model_data, function(x) {
  dat <- t(x[["vote"]])
  splicing_gic_logistic <- slide(dat, weight = rep(1, nrow(dat)), 
                                 tune.type = "gic", ic.scale = 1)
  splicing_gic_logistic[[1]]
})
save(senate_nodewise, file = "senate_nodewise.rda")

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
  g %v% "text" <- ""
  
  ggnetwork_default <- ggnetwork::ggnetwork(g, layout = "eigen")
  return(ggnetwork_default)
}
load("senate_nodewise.rda")
load("senate_model_data.rda")

congress_range <- 112:117
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
  facet_wrap(congress ~ ., scales = "free", ncol = 3) + 
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
pic
