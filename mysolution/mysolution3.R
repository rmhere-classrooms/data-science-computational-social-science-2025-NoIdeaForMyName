library(shiny)
library(igraph)
library(bslib)
library(memoise)

# Logika aplikacji

loadGraph <- function() {
  dfGraph <- read.csv2(
    url("http://bergplace.org/share/out.radoslaw_email_email"),
    skip = 2, sep = " "
  )[, 1:2]
  
  g <- graph.data.frame(dfGraph, directed = TRUE)
  g <- simplify(g)
  
  dfGraph <- dfGraph[dfGraph$X1 != dfGraph$X2, ]
  
  cat("Liczba węzłów w grafie:", vcount(g), "\n")
  cat("Liczba krawędzi w grafie:", ecount(g), "\n")
  
  dfSummary <- aggregate(
    rep(1, nrow(dfGraph)), # liczymy po jednym dla każdego wiersza
    by = list(X1 = dfGraph$X1, X2 = dfGraph$X2),
    FUN = sum
  )
  names(dfSummary)[3] <- "emailNb"
  
  # Sumy dla każdego wierzchołka
  vertex_totals <- aggregate(
    emailNb ~ X1, 
    data = dfSummary, 
    FUN = sum
  )
  names(vertex_totals)[2] <- "total_emails"

  email_lookup <- dfSummary
  vertex_lookup <- vertex_totals
  
  edges <- as_edgelist(g)
  
  # Wektor wag
  weights <- numeric(ecount(g))
  
  for (i in 1:ecount(g)) {
    from_node <- edges[i, 1]
    to_node <- edges[i, 2]
    
    # Liczba maili między określonymi węzłami
    email_count <- email_lookup[email_lookup$X1 == from_node & email_lookup$X2 == to_node, "emailNb"]

    # Liczba maili z węzła wychodzącego
    total_emails <- vertex_lookup[vertex_lookup$X1 == from_node, "total_emails"]

    weights[i] <- email_count / total_emails
  }
  
  E(g)$weight <- weights
  
  # Test poprawności przypisanych wag
  cat("Test poprawności...", "\n")
  for (v in V(g)) {
    v_name <- V(g)$name[v]
    outgoing <- incident(g, v, mode = "out")
    if (length(outgoing) > 0) {
      total_weight <- sum(E(g)$weight[outgoing])
      if (total_weight < 0.99999) { # Poprawka ze względu na błędy zmiennoprzecinkowe
        cat("TEST FAILED: vertex", v_name, "sum =", total_weight, "\n")
      }
    }
  }
  cat("Koniec testu...", "\n")
  
  return(g)
}

tryInfect <- function(w_ij, activateProbability) {
  return(w_ij * activateProbability > runif(1))
}

get_neighbors_list_helper <- (function (g, n) lapply(1:n, function(i) neighbors(g, i, mode = "out")))
get_neighbors_list <- memoise(get_neighbors_list_helper)

get_edge_weights_helper <- function(g, n) {
  edge_weights <- matrix(0, n, n)
  edges <- as_edgelist(g, names = FALSE)
  weights <- E(g)$weight
  for (i in seq_len(nrow(edges))) {
    edge_weights[edges[i,1], edges[i,2]] <- weights[i]
  }
  return(edge_weights)
}
get_edge_weights <- memoise(get_edge_weights_helper)

spreadIndependentCascades <- function(g, initialActivated, activateProbability = 1, iterationsNb = 10) {
  # activateProbability [0.1; 2]
  # activated - zarażony w ostatniej turze i może teraz zarażać
  # actived - zarażająy w tej turze - od kolejnej nie może już zarażać
  
  n <- vcount(g)
  
  activated <- logical(n)
  actived <- logical(n)
  activated[initialActivated] <- TRUE

  neighbors_list <- get_neighbors_list(g, n)
  edge_weights <- get_edge_weights(g, n)

  activatedList <- initialActivated
  activatedPerIteration <- integer(0)
  for (iter in 1:iterationsNb) {
    activatedNb <- 0
    oldActivatedList <- activatedList
    activatedList <- integer(0)
    for (spreader in oldActivatedList) {
      actived[spreader] <- T
      neighbors_out <- neighbors_list[[spreader]]
      for (n in neighbors_out) {
        if (!activated[n] && !actived[n]) {
          w_ij <- edge_weights[spreader, n]
          if (tryInfect(w_ij, activateProbability)) {
            activated[n] <- T
            activatedList <- c(activatedList, n)
            activatedNb <- activatedNb + 1
          }
        }
      }
    }

    activatedPerIteration[iter] <-  activatedNb
  }
  
  return(activatedPerIteration)
}

# Strategie wyboru węzłów początkowych
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

# V - reversed max pageRank
# Odwrócony PageRank oblicza klasyczną miarę, ale na grafie z odwróconymi 
# kierunkami krawędzi, dzięki czemu wskazuje węzły mające duży wpływ 
# „na zewnątrz”, czyli silnie oddziałujące na inne. W takim ujęciu wysoką 
# wartość otrzymują te węzły, które w oryginalnym grafie mają wiele lub ważne 
# połączenia wychodzące, co czyni je dobrymi starterami do rozprzestrzeniania 
# informacji. Dzięki temu miara lepiej pasuje do modeli dyfuzji w sieciach 
# kierunkowych niż klasyczny PageRank.
chooseReversedMaxPagerank <- function(g, nb) {
  pr <- page_rank(reverse_edges(g), directed = TRUE)$vector
  top_indices <- order(pr, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}

# EXPERIMENTING...
experimentSpreading <- function(g, initialNb, initialChooseFunction, n = 100, maxIter = 10, activProb) {
  infectedIterationsList = list()
  maxLength <- 0
  for (i in 1:n) {
    result <- spreadIndependentCascades(
      g = g, 
      initialActivated = initialChooseFunction(g, initialNb), 
      activateProbability = activProb, 
      iterationsNb = maxIter
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

runFullExperiment <- function(g, initialNb, activProb, iterNb, n = 100) {
  chooseFunctionList <- list(
    "max degree"          = chooseMaxDegree,
    "max betweenness"     = chooseMaxBetweenness,
    "max closeness"       = chooseMaxCloseness,
    "random"              = chooseRandom,
    "rev max pagerank"    = chooseReversedMaxPagerank
  )
  
  results <- list()
  for (name in names(chooseFunctionList)) {
    cat("Eksperyment dla", name, "\n")
    results[[name]] <- experimentSpreading(
      g = g,
      initialNb = initialNb,
      initialChooseFunction = chooseFunctionList[[name]],
      n = n,
      maxIter = iterNb,
      activProb = activProb
    )
  }
  
  return(results)
}

# ShinyApp User Interface

ui <- page_sidebar(
  title = "Rozprzestrzenianie się informacji w sieciach",
  sidebar = sidebar(
    sliderInput(
      inputId = "activProbId",
      label = "Prawdopodobieństwo aktywacji:",
      min = 10,
      max = 200,
      value = 100,
      post = "%"
    ),
    
    sliderInput(
      inputId = "iterNbId",
      label = "Liczba iteracji:",
      min = 1,
      max = 50,
      value = 10
    )
  ),
  plotOutput(outputId = "diffusionPlot")
)

# ShinyApp Server

server <- function(input, output, session) {
  
  g <- loadGraph()   # początkowe ładowanie grafu
  
  output$diffusionPlot <- renderPlot({
    
    activateProb <- input$activProbId / 100
    iterNb <- input$iterNbId
    initialNb <- max(1, round(0.05 * vcount(g)))  # 5% wierzchołków
    
    results <- runFullExperiment(
      g = g,
      initialNb = initialNb,
      activProb = activateProb,
      iterNb = iterNb,
      n = 100
    )
    
    allValues <- unlist(results)
    plot(
      1:iterNb, 
      results[[1]], 
      type = "l", 
      lwd = 2,
      ylim = c(min(allValues), max(allValues)),
      xlab = "Iteracja", ylab = "Aktywowane węzły",
      main = "Dyfuzja informacji"
    )
    
    colIndex <- 2
    for (name in names(results)[-1]) {
      lines(
        1:iterNb, 
        results[[name]], 
        lwd = 2, 
        col = colIndex
      )
      colIndex <- colIndex + 1
    }
    
    legend(
      "topright",
      legend = names(results),
      col = 1:length(results),
      lwd = 2
    )
  })
}

shinyApp(ui = ui, server = server)
