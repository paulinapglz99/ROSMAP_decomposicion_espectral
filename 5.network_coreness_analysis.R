#Coreness analysis

#Libraries ----

pacman::p_load('tidyverse', 
               'igraph')

#Read graph ----

graph <- read.graph(file = "noMCI_17112023_graph_MI0.532.graphml",
                    format = "graphml")

#coreness ----

coreness <- coreness(graph) %>% 
  as.data.frame() %>% 
  rename("core_by_node" = ".") 

coreness$gene <- rownames(coreness)

rownames(coreness) <- NULL

#Histogram of the coreness distribution

hist(coreness$core_by_node)

max(coreness$core_by_node) #el maximo de coreness por nodo para este grafo

hist_coreness <- ggplot(coreness, aes(x = core_by_node)) +
  geom_histogram(fill = "skyblue", color = "white") +
  labs(title = "Histograma de coreness por nodo",
       subtitle = "Para sujetos ROSMAP sin demencia,despues del corte de percolacion de MI",
       x = "Coreness by node", y = "Frecuencia") +
  scale_x_continuous(breaks = seq(0, 160, by = 10)) +
  theme_light()

# Filtering by coreness -----

coreness_filter <- coreness %>% 
  filter(core_by_node >= 20) #filtrar para tener un grafo con nodos con un core max de __

hist(coreness_filter$core_by_node)

#Histogram of filtered coreness----

hist_coreness_filter <- ggplot(coreness_filter, aes(x = core_by_node)) +
  geom_histogram(fill = "skyblue", color = "white") +
  labs(title = "Histograma de coreness por nodo con filtro",
       subtitle = "Para sujetos ROSMAP sin demencia, despues del corte de percolacion de MI",
       x = "Coreness by node", y = "Frecuencia") +
  scale_x_continuous(breaks = seq(0, 160, by = 10)) +
  theme_light()

#Now the subgraph ----

# Create a vector with the name of nodes I want to subset

coreness_filter_v <- as.vector(coreness_filter$gene)

coreness_filter_v <- unlist(coreness_filter_v)

#haciendo un subgrafo solo con los vertices del core we want

subgraph_kcore <- induced_subgraph(graph,
                                    vids = coreness_filter_v)

#plotear 

subraph_kcore_filter.p <- plot(subgraph_kcore)
