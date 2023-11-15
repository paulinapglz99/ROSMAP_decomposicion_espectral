#Analisis de k-cores para redes

#paquetes

pacman::p_load('tidyverse', 
               'igraph',
               'vroom')

#llamar a la red

matrix <- vroom(file = '/datos/rosmap/matriz_coexpre_20231011.txt')

#subset pa aprender ---

matrix_subset <- matrix[1:100, 1:100] %>% 
  pivot_longer(cols = -gene, names_to = 'gene_to', values_to = 'mi') %>% 
  filter(mi >= 0.5) 

#hacer red a partir de data frame

graph <- graph_from_data_frame(matrix_subset, 
                      directed = F)  #por fin veo una red :'u

#wacharla --

plot(graph)

#generalidades de la red ---

components(graph) #saber cuantos componentes tiene

degree(graph) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

#coreness ----

coreness <- coreness(graph) %>% 
  as.data.frame()

coreness_3 <- coreness %>% 
  filter(. == 3) #filtrar para tener solo esos cores

coreness_3 <- as.vector(rownames(coreness_3))

#haciendo un subgrafo solo con los vertices de coreness_3
subraph_kcore3 <- induced_subgraph(graph, vids = coreness_3)

#plotear 

print(subraph_kcore3)
subraph_kcore3.p <- plot(subraph_kcore3) #IFUCKINGDIDIT
save.image()


###sin plotear para que no muera esta cosa
#hacer red a partir de data frame

graph_complete <- graph_from_data_frame(matrix, 
                               directed = F)  #por fin veo una red :'u


#generalidades de la red ---

components(graph_complete) #saber cuantos componentes tiene

degree(graph_complete) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

#coreness ----

coreness <- coreness(graph_complete) %>% 
  as.data.frame()

coreness_3 <- coreness %>% 
  filter(. == 3) #filtrar para tener solo esos cores

coreness_3 <- as.vector(rownames(coreness_3))

#haciendo un subgrafo solo con los vertices de coreness_3
subraph_kcore3 <- induced_subgraph(graph, vids = coreness_3)



