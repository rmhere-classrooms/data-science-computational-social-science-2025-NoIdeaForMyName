"
1. Plik z rozwiązaniem nazwij mysolution/mysolution2.R (lub inne rozszerzenie, zależnie
od języka programowania).
2. Wygeneruj graf wedle modelu Barabási-Albert z tysiącem węzłów
3. Zwizualizuj graf layoutem Fruchterman & Reingold
4. Znajdź najbardziej centralny węzeł według miary betweenness, jaki ma numer?
5. Jaka jest średnica grafu?
6. W komentarzu napisz czym różnią się grafy Barabási-Albert i Erdős-Rényi.
"

library(igraph)

# 2.
g <- barabasi.game(1000)

# 3.
layout <- layout.fruchterman.reingold(g)
plot(g, layout=layout, vertex.size=2,
     vertex.label=NA, edge.arrow.size=.2)

# 4.
most_central_v <- V(g)[betweenness(g)==max(betweenness(g))]
cat("Najbadziej centralny węzeł według miary betweenness ma numer:", most_central_v, "\n")

# 5.
cat("Średnica grafu to:", diameter(g))

# 6.
"
Grafy Erdős-Rényi i Barabási-Albert różnią się zasadniczo sposobem generowania 
połączeń między węzłami. 

W modelu Erdős-Rényi każda możliwa krawędź powstaje losowo i niezależnie 
z jednakowym prawdopodobieństwem, co prowadzi do grafów o równomiernym 
rozkładzie stopni węzłów i niskiej klasteryzacji. 

Z kolei grafy Barabási-Albert powstają w procesie preferencyjnego przyłączania 
nowych węzłów do już dobrze połączonych węzłów (takich o wysokim stopniu), 
co skutkuje powstaniem tzw. hubów.

Podsumowując, w sieciach Barabási-Albert występują centralne, 
silnie połączone węzły, podczas gdy grafy Erdős-Rényi są bardziej jednorodne 
i nie mają wyraźnych hubów.
"
