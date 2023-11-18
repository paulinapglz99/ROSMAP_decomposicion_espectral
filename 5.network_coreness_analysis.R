#Coreness analysis

#Libraries ----

pacman::p_load('tidyverse', 
               'igraph')

#Read graph ----

graph <- read.graph(file = "noMCI_17112023_graph_MI0.532.graphml",
                    format = "graphml")

#coreness ----

coreness <- coreness(graph) %>% 
  as.data.frame()

#Histogram of the coreness distribution

hist(coreness$.)

# Filtering by coreness -----

coreness_25 <- coreness %>% 
  filter(. >= 25) #filtrar para tener un grafo con nodos con un core max de 25

hist(coreness_25$.)

#Now the subgraph ----

# Create a vector with the name of nodes I want to subset

coreness_25v <- as.vector(rownames(coreness_25))

#haciendo un subgrafo solo con los vertices de coreness_3

subraph_kcore3 <- induced_subgraph(graph, vids = coreness_25)

#plotear 

print(subraph_kcore3)
subraph_kcore3.p <- plot(subraph_kcore3) #IFUCKINGDIDIT
save.image()


###sin plotear para que no muera esta cosa
#hacer red a partir de data frame

graph_complete <- graph_from_data_frame(matrix, 
                               directed = F)  #por fin veo una red :'u
#coreness ----

coreness <- coreness(graph_complete) %>% 
  as.data.frame()

coreness_3 <- coreness %>% 
  filter(. == 3) #filtrar para tener solo esos cores

coreness_3 <- as.vector(rownames(coreness_3))

#haciendo un subgrafo solo con los vertices de coreness_3

subraph_kcore3 <- induced_subgraph(graph, vids = coreness_3)
