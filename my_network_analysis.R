#analisis de redes

#paquetes

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#llamar a la matriz de adyacencia

matrix <- vroom::vroom(file = '/datos/rosmap/matriz_coexpre_allAD_11052023_zero.txt') %>% 
   as.data.frame()

#antes de hacer la red ----

matrix_MI <- matrix %>% 
  pivot_longer(cols = -gene,
               names_to = 'gene_to',
               values_to = 'mi')

########---

#histograma de Mutual Information


hist(matrix_MI$mi)

hist_allAD_MI <- ggplot(matrix_MI, aes(x = mi)) +
  geom_histogram() +
  labs(title = "Histograma de Informacion Mutua",
       subtitle = "Para sujetos ROSMAP con cogdx de AD",
       x = "Informacion Mutua", y = "Frecuencia") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))

#decision, voy a cortar en 0.5 y en 1 para ver como queda

matrix_MI_0.5 <- matrix_MI %>% 
 filter(mi >= 0.5)

hist(matrix_MI_0.5$mi)


matrix_MI_1.0 <- matrix_MI %>% 
  filter(mi >= 1)

hist(matrix_MI_1.0$mi)


#100 000 edges de mayor a menor y me quedo con los primeros 10 000 interacciones
#quedate con el 5% de los nodos mas altos

#hacer red a partir de data frame

graphAD <- graph_from_data_frame(matrix_MI, 
                                 directed = F)

#graphAD_0.5 <- graph_from_data_frame(matrix_MI_0.5, 
#                               directed = F)

#graphAD_1.0 <- graph_from_data_frame(matrix_MI_1.0, 
#                                   directed = F)  


write_graph(graphAD, 
           file = 'coexpre_allAD_11052023_graph_nocutt.txt', 
           'edgelist')

vroom::vroom(file = "coexpre_allAD_11052023_graph_AD.txt")

#generalidades de las redes ---

components(graphAD)

components(graph_0.5) #saber cuantos componentes tiene

components(graph_1.0) #saber cuantos componentes tiene
#we can also see degree by node
degree(graphAD)

degree(graph_0.5) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

degree(graph_1.0) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

degree_distribution(graphAD) %>%
  unique()

#edge betweenneess -----

betweenneessAD <- betweenness(graphAD,
            directed = F)

#coreness ----

corenessAD <- coreness(graphAD) %>%  
  as.data.frame()

coreness_0.5 <- coreness(graph_0.5) %>%  
  as.data.frame()

coreness_1.0 <- coreness(graph_1.0) %>%  
  as.data.frame()

hist(corenessAD$.)
hist(coreness_0.5$.)  #histograma de mi coreness distribucion de probabilidad

hist(coreness_1.0$.)  #histograma de mi coreness distribucion de probabilidad

#extraigo vectores

coreness_0.5v <- as.vector(rownames(coreness_0.5))

coreness_1.0v <- as.vector(rownames(coreness_1.0))

#haciendo un subgrafo solo con los vertices de coreness_3

subraph_kcore <- induced_subgraph(graph, vids = coreness_1.0)

#plotear 

print(subraph_kcore3)
subraph_kcore3.p <- plot(subraph_kcore3) #IFUCKINGDIDIT
