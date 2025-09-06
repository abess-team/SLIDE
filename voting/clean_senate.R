rm(list = ls()); gc(reset = TRUE)
library(reshape2)

path <- "~/voting"
setwd(path)

congress_study <- 106:117
year <- seq(2001, 2023, by = 2)

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
# votes <- votes[votes[["cast_code"]] %in% c(1, 6, 9), ]
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
