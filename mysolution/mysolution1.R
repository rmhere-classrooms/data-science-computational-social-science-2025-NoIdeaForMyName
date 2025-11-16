# 1. Plik z rozwiązaniem nazwij mysolution/mysolution1.R (lub inne rozszerzenie, zależnie
# od języka programowania).
# 2. Wygeneruj sieć Erdős-Rényi o stu wierzchołkach i prawdopodobieństwie krawędzi = 0.05.
# 3. Wydrukuj podsumowanie grafu - czy graf jest ważony?
# 4. Wylistuj wszystkie wierzchołki i krawędzie.
# 5. Ustaw wagi wszystkich krawędzi na losowe z zakresu 0.01 do 1
# 6. Wydrukuj ponownie podsumowanie grafu - czy teraz graf jest ważony?
# 7. Jaki jest stopień każdego węzła? Następnie stwórz histogram stopni węzłów.
# 8. Ile jest klastrów (connected components) w grafie?
# 9. Zwizualizuj graf w taki sposób, aby rozmiar węzłów odpowiadał mierze PageRank.


library(igraph)

# 2.
g <- erdos.renyi.game(p.or.m=0.05, n=100)

# 3.
summary(g)

# Graf nie jest ważony - podsumowanie grafu to:

# IGRAPH c43c49b U--- 100 241 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n)

# Gdyby graf był ważony, na poz. 20. w 1. wierszu byłoby 'W' zamiast '-'.


# 4.
V(g)
E(g)

# 5.
E(g)$weight <- runif(length(E(g)), 0.01, 1)

# 6.
summary(g)
# Po tej operacji graf jest ważony:

# IGRAPH c43c49b U-W- 100 213 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n), weight (e/n)


# 7.
degree(g)
hist(degree(g))

# 8.
cl <- clusters(g)
cat("Liczba klastrów w grafie:", cl$no, "\n")

# 9.
pr <- page.rank(g)$vector
plot(g, vertex.size=pr*300, vertex.label=NA, edge.arrow.size=.2)
