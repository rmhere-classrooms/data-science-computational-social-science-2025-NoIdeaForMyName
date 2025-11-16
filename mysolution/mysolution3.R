library(shiny)
library(igraph)
library(bslib)
library(memoise)


# Logika aplikacji

loadGraph <- function() {
  df_graph <- read.csv2(
    url("http://bergplace.org/share/out.radoslaw_email_email"),
    skip = 2, sep = " "
  )[, 1:2]
  
  g <- graph.data.frame(df_graph, directed = TRUE)
  g <- simplify(g)
  
  df_graph <- df_graph[df_graph$X1 != df_graph$X2, ]
  
  #cat("Liczba węzłów w grafie:", vcount(g), "\n")
  #cat("Liczba krawędzi w grafie:", ecount(g), "\n")
  
  df_summary <- aggregate(
    rep(1, nrow(df_graph)), # liczymy po jednym dla każdego wiersza
    by = list(X1 = df_graph$X1, X2 = df_graph$X2),
    FUN = sum
  )
  names(df_summary)[3] <- "email_nb"
  
  # Sumy dla każdego wierzchołka
  vertex_totals <- aggregate(
    email_nb ~ X1, 
    data = df_summary, 
    FUN = sum
  )
  names(vertex_totals)[2] <- "total_emails"

  email_lookup <- df_summary
  vertex_lookup <- vertex_totals
  
  edges <- as_edgelist(g)
  
  # Wektor wag
  weights <- numeric(ecount(g))
  
  for (i in 1:ecount(g)) {
    from_node <- edges[i, 1]
    to_node <- edges[i, 2]
    
    # Liczba maili między określonymi węzłami
    email_count <- email_lookup[email_lookup$X1 == from_node & email_lookup$X2 == to_node, "email_nb"]

    # Liczba maili z węzła wychodzącego
    total_emails <- vertex_lookup[vertex_lookup$X1 == from_node, "total_emails"]

    weights[i] <- email_count / total_emails
  }
  
  E(g)$weight <- weights
  
  # Test poprawności przypisanych wag
  #cat("Test poprawności...", "\n")
  for (v in V(g)) {
    v_name <- V(g)$name[v]
    outgoing <- incident(g, v, mode = "out")
    if (length(outgoing) > 0) {
      total_weight <- sum(E(g)$weight[outgoing])
      if (total_weight < 0.99999) { # Poprawka ze względu na błędy zmiennoprzecinkowe
        warning("WEIGHT TEST FAILED: vertex", v_name, "sum =", total_weight, "\n")
      }
    }
  }
  #cat("Koniec testu...", "\n")
  
  return(g)
}

tryInfect <- function(w_ij, activate_probability) {
  return(w_ij * activate_probability > runif(1))
}

getNeighborsListHelper <- (function (g, n) lapply(1:n, function(i) neighbors(g, i, mode = "out")))
getNeighborsList <- memoise(getNeighborsListHelper)

getEdgeWeightsHelper <- function(g, n) {
  edge_weights <- matrix(0, n, n)
  edges <- as_edgelist(g, names = FALSE)
  weights <- E(g)$weight
  for (i in seq_len(nrow(edges))) {
    edge_weights[edges[i,1], edges[i,2]] <- weights[i]
  }
  return(edge_weights)
}
getEdgeWeights <- memoise(getEdgeWeightsHelper)

spreadIndependentCascades <- function(g, initial_activated, activate_probability = 1, iterations_nb = 10) {
  # activate_probability [0.1; 2]
  # activated - zarażony w ostatniej turze i może teraz zarażać
  # actived - zarażająy w tej turze - od kolejnej nie może już zarażać
  
  n <- vcount(g)
  
  activated <- logical(n)
  actived <- logical(n)
  activated[initial_activated] <- TRUE

  neighbors_list <- getNeighborsList(g, n)
  edge_weights <- getEdgeWeights(g, n)

  activated_list <- initial_activated
  activated_per_iteration <- integer(0)
  for (iter in 1:iterations_nb) {
    activated_nb <- 0
    oldactivated_list <- activated_list
    activated_list <- integer(0)
    for (spreader in oldactivated_list) {
      actived[spreader] <- T
      neighbors_out <- neighbors_list[[spreader]]
      for (n in neighbors_out) {
        if (!activated[n] && !actived[n]) {
          w_ij <- edge_weights[spreader, n]
          if (tryInfect(w_ij, activate_probability)) {
            activated[n] <- T
            activated_list <- c(activated_list, n)
            activated_nb <- activated_nb + 1
          }
        }
      }
    }

    activated_per_iteration[iter] <-  activated_nb
  }
  
  return(activated_per_iteration)
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
# kierunkami krawędzi, dzięki czemu wskazuje na węzły mające wiele połączeń
# na zewnątrz (do innych "ważnych" węzłów, które w odwróconej wersji 
# wskazywały na te węzły). Ta właściwość czyni je dobrymi punktami 
# początkowymi do rozprzestrzeniania informacji. Dzięki temu miara lepiej 
# pasuje do modeli dyfuzji w sieciach kierunkowych niż klasyczny PageRank
chooseReversedMaxPagerank <- function(g, nb) {
  pr <- page_rank(reverse_edges(g), directed = TRUE)$vector
  top_indices <- order(pr, decreasing = TRUE)[1:nb]
  return(V(g)[top_indices])
}

# EXPERIMENTING...
experimentSpreading <- function(g, initial_nb, initialChooseFunction, n = 100, max_iter = 10, activ_prob) {
  infected_iterations_list = list()
  max_length <- 0
  for (i in 1:n) {
    result <- spreadIndependentCascades(
      g = g, 
      initial_activated = initialChooseFunction(g, initial_nb), 
      activate_probability = activ_prob, 
      iterations_nb = max_iter
    )
    max_length <- max(max_length, length(result))
    while (length(result) < max_iter) {
      result <- c(result, 0)
    }
    infected_iterations_list[[i]] <- result
  }
  summed <- Reduce(`+`, infected_iterations_list)
  final_result <- summed / n
  return(final_result)
}

runFullExperiment <- function(g, initial_nb, activ_prob, iter_nb, n = 100) {
  chooseFunctionList <- list(
    "max degree"          = chooseMaxDegree,
    "max betweenness"     = chooseMaxBetweenness,
    "max closeness"       = chooseMaxCloseness,
    "random"              = chooseRandom,
    "rev max pagerank"    = chooseReversedMaxPagerank
  )
  
  results <- list()
  for (name in names(chooseFunctionList)) {
    #cat("Eksperyment dla", name, "\n")
    results[[name]] <- experimentSpreading(
      g = g,
      initial_nb = initial_nb,
      initialChooseFunction = chooseFunctionList[[name]],
      n = n,
      max_iter = iter_nb,
      activ_prob = activ_prob
    )
  }
  
  return(results)
}

# ShinyApp User Interface

ui <- page_sidebar(
  title = "Rozprzestrzenianie się informacji w sieciach",
  sidebar = sidebar(
    sliderInput(
      inputId = "activ_probId",
      label = "Prawdopodobieństwo aktywacji:",
      min = 10,
      max = 200,
      value = 100,
      post = "%"
    ),
    
    sliderInput(
      inputId = "iter_nbId",
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
    
    activate_prob <- input$activ_probId / 100
    iter_nb <- input$iter_nbId
    initial_nb <- max(1, round(0.05 * vcount(g)))  # 5% wierzchołków
    
    results <- runFullExperiment(
      g = g,
      initial_nb = initial_nb,
      activ_prob = activate_prob,
      iter_nb = iter_nb,
      n = 100
    )
    
    all_values <- unlist(results)
    plot(
      1:iter_nb, 
      results[[1]], 
      type = "l", 
      lwd = 2,
      col = 1,
      ylim = c(min(all_values), max(all_values)),
      xlab = "Iteracja", ylab = "Aktywowane węzły",
      main = "Dyfuzja informacji"
    )
    points(
      1:iter_nb, 
      results[[1]], 
      pch = 19,
      col = 1
    )
    
    col_index <- 2
    for (name in names(results)[-1]) {
      lines(
        1:iter_nb, 
        results[[name]], 
        lwd = 2, 
        col = col_index
      )
      points(
        1:iter_nb,
        results[[name]],
        col = col_index,
        pch = 19
      )
      col_index <- col_index + 1
    }
    
    legend(
      "topright",
      legend = names(results),
      col = 1:length(results),
      lwd = 2,
      pch = 19
    )
  })
}

shinyApp(ui = ui, server = server)
