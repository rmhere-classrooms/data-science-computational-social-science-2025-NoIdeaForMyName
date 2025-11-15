library(igraph)


# 4.
setwd("/Users/mmaksanty/Documents/Studia/7 semestr/danologia/lista_5/solution/data-science-computational-social-science-2025-NoIdeaForMyName")
dfGraph <- read.csv2("css/out.radoslaw_email_email", skip=2, sep= " ")[, 1:2]
g <- graph.data.frame(dfGraph, directed = T)


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


for (i in seq_along(V(g))) {
  v <- V(g)[i]
  cnt_i <- sum(dfSummary$emailNb[dfSummary$X1 == v$name])
  test <- 0
  neigh <- neighbors(g, v, mode = "out")
  for (j in seq_along(neigh)) {
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


# 7.
tryInfect <- function(g, v1, v2, activateProbability) {
  w_ij <- E(g)[v1 %->% v2]$weight
  return(w_ij * activateProbability > runif(1))
}

spreadIndependentCascades <- function(g, initialActivated, activateProbability = 1, maxIterations = 10) {
  # activateProbability [0.1; 2]
  # activated - zarażony w ostatniej turze i może teraz zarażać
  # actived - zarażająy w tej turze - od kolejnej nie może już zarażać
  V(g)$activated <- F
  V(g)$actived <- F
  V(g)[initialActivated]$activated <- T
  
  activatedList <- V(g)[initialActivated]
  newlyActivated <- T
  iterationsLeft <- maxIterations
  activatedPerIteration <- c()
  while(iterationsLeft > 0 && newlyActivated) {
    #cat("Iteration:", maxIterations-iterationsLeft+1, "\n")
    activatedNb <- 0
    newlyActivated <- F
    oldActivatedList <- activatedList
    activatedList <- V(g)[c()]
    #cat("ACTIVATION LIST:", oldActivatedList, "\n")
    for (spreader in oldActivatedList) {
      V(g)[spreader]$actived <- T
      neighbors_out <- neighbors(g, spreader, mode = "out")
      for (n in neighbors_out) {
        if (!V(g)[n]$activated && !V(g)[n]$actived && tryInfect(g, spreader, n, activateProbability)) {
          #cat("spreading!", "\n")
          V(g)[n]$activated <- T
          activatedList <- append(activatedList, n)
          activatedNb <- activatedNb + 1
          newlyActivated <- T
        }
      }
    }
    
    iterationsLeft <- iterationsLeft - 1
    activatedPerIteration <- c(activatedPerIteration, activatedNb)
  }
  
  return(activatedPerIteration)
}


# 8.

# I - max outdegree
chooseMaxDegree <- function(g, nb) {
  outdeg <- degree(g, mode = "out")
  top_indices <- order(outdeg, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}

# II - max betweeness
chooseMaxBetweenness <- function(g, nb) {
  btw <- betweenness(g)
  top_indices <- order(btw, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}

# III - max closeness
chooseMaxCloseness <- function(g, nb) {
  cls <- closeness(g, mode = "out")
  top_indices <- order(cls, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}

# IV - random vertices
chooseRandom <- function(g, nb) {
  return(sample(V(g), nb))
}

# V - max pageRank
# TODO opis miary pageRank
chooseMaxPagerank <- function(g, nb) {
  pr <- page_rank(g, directed = TRUE)$vector
  top_indices <- order(pr, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}


# EXPERIMENTING...
experimentSpreading <- function(g, initialNb, initialChooseFunction, n = 100, maxIter = 10) {
  infectedIterationsList = list()
  maxLength <- 0
  for (i in 1:n) {
    result <- spreadIndependentCascades(
      g = g, 
      initialActivated = initialChooseFunction(g, initialNb), 
      activateProbability = 1, 
      maxIterations = maxIter
    )
    maxLength <- max(maxLength, length(result))
    while (length(result) < maxIter) {
      result <- c(result, 0)
    }
    infectedIterationsList[[i]] <- result
  }
  summed <- Reduce(`+`, infectedIterationsList)
  finalResult <- summed / n
  return(finalResult)
}

fullExperiment <- function(g, initialNb, chooseFunctionList, n = 5, maxIter = 10) {
  
  results <- list()
  
  for (name in names(chooseFunctionList)) {
    cat("Eksperyment dla", name, "\n")
    chooseFunc <- chooseFunctionList[[name]]
    result <- experimentSpreading(
      g = g,
      initialNb = initialNb,
      initialChooseFunction = chooseFunc,
      n = n,
      maxIter = maxIter
    )
    results[[name]] <- result
  }
  
  # Rysowanie wykresu (1. seria danych)
  allValues <- unlist(results)
  minY <- min(allValues)
  maxY <- max(allValues)
  
  firstName <- names(results)[1]
  firstData <- results[[firstName]]
  
  plot(
    x = 1:length(firstData),
    y = firstData,
    type = "l",
    col = 1,
    lwd = 2,
    xlab = "Iteracja",
    ylab = "Liczba aktywowanych węzłów",
    main = "Przebieg dyfuzji informacji",
    ylim = c(minY, maxY) 
  )
  
  # Kolejne serie
  i <- 2
  for (name in names(results)[-1]) {
    lines(
      1:length(results[[name]]), 
      results[[name]], 
      col = i, 
      lwd = 2
    )
    i <- i + 1
  }
  
  # Legenda
  legend(
    "topright",
    legend = names(results),
    col = 1:length(results),
    lwd = 2,
  )
  
}

initialInfectedNb <- 0.05 * vcount(g) # 5% wszystkich wierzchołków
chooseFunctionList <- list(
  "max degree" = chooseMaxDegree,
  "max betweenness" = chooseMaxBetweenness,
  "max closeness" = chooseMaxCloseness,
  "random" = chooseRandom,
  "max pagerank" = chooseMaxPagerank
)

fullExperiment(
  g = g, 
  initialNb = initialInfectedNb, 
  chooseFunctionList = chooseFunctionList, 
  n = 5, 
  maxIter = 10
)

