library(shiny)
library(igraph)
library(bslib)

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
  return(g)
}

tryInfect <- function(g, v1, v2, activateProbability) {
  w_ij <- E(g)[v1 %->% v2]$weight
  return(w_ij * activateProbability > runif(1))
}

spreadIndependentCascades <- function(g, initialActivated, activateProbability = 1, iterationsNb = 10) {
  # activateProbability [0.1; 2]
  # activated - zarażony w ostatniej turze i może teraz zarażać
  # actived - zarażająy w tej turze - od kolejnej nie może już zarażać
  V(g)$activated <- F
  V(g)$actived <- F
  V(g)[initialActivated]$activated <- T
  
  activatedList <- V(g)[initialActivated]
  iterationsLeft <- iterationsNb
  activatedPerIteration <- c()
  while(iterationsLeft > 0) {
    #cat("Iteration:", iterationsNb-iterationsLeft+1, "\n")
    activatedNb <- 0
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
        }
      }
    }
    
    iterationsLeft <- iterationsLeft - 1
    activatedPerIteration <- c(activatedPerIteration, activatedNb)
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

runFullExperiment <- function(g, initialNb, activProb, iterNb, n = 5) {
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
      n = 1 # 5
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
