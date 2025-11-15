library(igraph)


# 4.
setwd("/Users/mmaksanty/Documents/Studia/7 semestr/danologia/lista_5/solution/data-science-computational-social-science-2025-NoIdeaForMyName")
dfGraph <- read.csv2("css/out.radoslaw_email_email", skip=2, sep= " ")[, 1:2]
g <- graph.data.frame(dfGraph, directed = TRUE)


# 5.
g <- simplify(g)
dfGraph <- dfGraph[dfGraph$X1 != dfGraph$X2, ] # usunięcie pętli z grafu
cat("Liczba węzłów w grafie:", vcount(g))
cat("Liczba krawędzi w grafie:", ecount(g))


# 6.
dfSummary <- aggregate(
  rep(1, nrow(dfGraph)), # liczymy po jednym dla każdego wiersza
  by = list(X1 = dfGraph$X1, X2 = dfGraph$X2),
  FUN = sum
)
names(dfSummary)[3] <- "emailNb"
head(dfSummary)


for(i in seq_along(V(g))) {
  v <- V(g)[i]
  cnt_i <- sum(dfSummary$emailNb[dfSummary$X1 == v$name])
  test <- 0
  neigh <- neighbors(g, v, mode = "out")
  for(j in seq_along(neigh)) {
    v_n <- neigh[j]
    cnt_ij <- dfSummary$emailNb[dfSummary$X1 == v$name & dfSummary$X2 == v_n$name]
    w_ij <- cnt_ij / cnt_i
    E(g)[v %->% v_n]$weight <- w_ij
    test <- test + w_ij
  }
  if (test < 0.99 && test > 0.01) { # test poprawności obliczeń
    cat("TEST:", test, "; ")
    cat("V I", i, "V", v, "N", neigh, "\n")
  }
}

