#
#7.topological_comparison.R
#This script makes the topological comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---
#install.packages("svglite")
pacman::p_load("igraph", 
               "ggraph",
               "tidyverse", 
               "gridExtra", 
               "svglite", 
               "tidygraph")

#Set seed --- --- 

set.seed(10)

#Declare functions --- --- 

#Declare function that compares edges

jaccard_edges <- function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

#Declare function that compares nodes

# #jaccard_nodes <- function(g1,g2){
#   a = sort(vertex.attributes(graph = g1)[["name"]])
#   b = sort(vertex.attributes(graph = g2)[["name"]])
#   
#   deLaCueva = length(intersect(a,b))/length(union(a,b))
#   return(deLaCueva)
# }

#Distribution 

process_distribution <- function(degree_distribution) {
  # Convert to data frame and count frequencies
  distribution.df <- as.data.frame(table(degree_distribution$degree))
  # Rename the column
  distribution.df <- distribution.df %>% rename(degree = "Var1")
  # Order by degree
  distribution.df <- distribution.df[order(distribution.df$degree, decreasing = F), ]
  # Calculate the cumulative sum
  distribution.df$CumulativeDegree <- cumsum(distribution.df$Freq)
  # Logarithm of the cumulative sum
  distribution.df$logCumulativeDegree <- log10(distribution.df$CumulativeDegree)
  # Add threshold
  CumulativeDegree_threshold <- quantile(distribution.df$logCumulativeDegree, probs = 0.95)
  # Puedes devolver el umbral si es necesario, por ejemplo:
  return(list(degree_distribution = distribution.df, cumm_threshold = CumulativeDegree_threshold))
}


#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                  graphnoAD = graphnoAD)

#Network topological comparison --- --- 

#Table of degree distribution

#Degree distributions of all graphs in my list
degree_distributions <- sapply(X = graphLists, FUN = degree)

#Build degree distribution dataframe

ADdegree <- degree_distributions[["graphAD"]]
ADdegree.df <- data.frame(gene = names(ADdegree), degree = ADdegree)
ADdegree_freq <- table(ADdegree.df$degree) %>% as.data.frame()
colnames(ADdegree_freq) <- c("degree", "Freq")
ADdegree_freq$degree <- as.numeric(as.character(ADdegree_freq$degree))
ADdegree_freq$Prob <- ADdegree_freq$Freq / sum(ADdegree_freq$Freq) # Frecuencia relativa

noADdegree <- degree_distributions[["graphnoAD"]]
noADdegree.df <- data.frame(gene = names(noADdegree), degree = noADdegree)
noADdegree_freq <- table(noADdegree.df$degree) %>% as.data.frame()
colnames(noADdegree_freq) <- c("degree", "Freq")
noADdegree_freq$degree <- as.numeric(as.character(noADdegree_freq$degree))
noADdegree_freq$Prob <- noADdegree_freq$Freq / sum(noADdegree_freq$Freq) # Frecuencia relativa

# Find the max value of "Freq" in both networks
max_freq_AD <- max(ADdegree_freq$Freq)
max_freq_noAD <- max(noADdegree_freq$Freq)

# Find max degree value in both graphs 
max_degree_AD <- max(as.numeric(as.character(ADdegree_freq$degree)))
max_degree_noAD <- max(as.numeric(as.character(noADdegree_freq$degree)))

# Define limits for make histograms comparable
max_y <- max(max_freq_AD, max_freq_noAD)
max_x <- max(max_degree_AD, max_degree_noAD)

#Plot both degree distributions

degree_disAD <- ggplot(ADdegree.df, aes(x = degree)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#961D4E", position = "identity") +
  labs(title = "Node degree distributions",
       subtitle = "for AD coexpression network", 
       x = "Degree",
       y = "Freq") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, max_x), breaks = seq(0, max_x, by = 10)) + 
  scale_fill_manual() +
  guides(fill = guide_legend(title = "Diagnosis"))

degree_disnoAD <- ggplot(noADdegree.df, aes(x = degree)) +
  geom_histogram(binwidth = 1,color = "black", fill = "#6153CC",  position = "identity") +
  labs(title = "Node degree distributions",
       subtitle = "for no AD coexpression network", 
       x = "Degree",
       y = "Freq") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +  # Normalizar el eje y
  scale_x_continuous(limits = c(0, max_x),breaks = seq(0, max_x, by = 10)) + 
  scale_fill_manual() +
  guides(fill = guide_legend(title = "Diagnosis"))

#Arrange in a grid

degree_dis <- grid.arrange(degree_disAD, degree_disnoAD, ncol =1 )

#Save plot
# 
# ggsave(filename = "bothdx_degree_distributions_coexpression_NIAReagan_histogram.jpg",
#       plot = degree_dis,
#       width = 25,
#       height = 20,
#       units = "cm",
#       dpi = 300,
#       )

#Calculate diameter of both graphs --- --- 

diameter <- sapply(X = graphLists, FUN = diameter)

# graphAD graphnoAD 
# 13       12 

#Eigenvector centrality of the network --- ---

eigen <- sapply(graphLists, FUN = eigen_centrality)

#Calculate clustering coefficient of both networks ---- ---

clustering_coefficient <- sapply(X = graphLists, FUN = transitivity)

# graphAD graphnoAD 
# 0.3357295 0.2979033 

infomap_modularity <- sapply(X = graphLists, FUN = cluster_infomap)

# Extract information from the modules

membership_modularity <- sapply(X = infomap_modularity, FUN = membership)

# Assign the modules as attributes of the vertices

V(graphAD)$modules <- membership_modularity[[1]]

V(graphnoAD)$modules <- membership_modularity[[2]]

#Comparison of modular structures between networks --- --- 

#To compare two networks at the modular level, it would be optimal to keep the same set of nodes. 
#Because the heuristic cut was made on the edges, we have an unequal number of nodes in the smaller
#network, so we will add the missing nodes to the smaller network. 

#graphs with equal node set

ADnodes <- V(graphAD)

noADnodes <- V(graphnoAD)

# Find elements in noADnodes but not in ADnodes
missing_elements_in_ADnodes <- setdiff(names(noADnodes), names(ADnodes))
length(missing_elements_in_ADnodes)

# Find elements in ADnodes but not in noADnodes
missing_elements_in_noADnodes <- setdiff(names(ADnodes), names(noADnodes))
length(missing_elements_in_noADnodes)

#Add nodes missing to AD graph

graphAD_plus <- add_vertices(graphAD, nv = length(missing_elements_in_ADnodes))

#Add nodes missing to noAD graph

graphnoAD_plus <- add_vertices(graphnoAD, nv = length(missing_elements_in_noADnodes))

#Set list of graphs

graphLists_plus <- list(graphAD_plus = graphAD_plus, 
                        graphnoAD_plus =graphnoAD_plus)


#Do they have similar nodes? --- --- 

#To make the topological comparison for nodes and edges we use the Jaccard index

#Apply function to compare nodes

NodesJaccard <- sapply(X = graphLists_plus, FUN = jaccard_nodes, g1 = graphLists_plus[[1]])
NodesJaccard

#graphAD graphnoAD  <-----
#0.6465551 1.0000000 

#Do they have similar edges? --- ---

#Apply function to compare edges

EdgesJaccard <- sapply(X = graphLists, FUN = jaccard_edges, g1 = graphLists[[1]])
EdgesJaccard

#Apply modularity algorithm --- ---

#Mutual information must be taken into account as weights of the edges

infomap_modularity <- list()

for (i in 1:length(graphLists)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  infomap_modularity[[i]] <- cluster_infomap(graph = graphLists[[i]], e.weights = graphLists[[i]]$mut_info_norm)
}

# Extract information from the modules

graphAD_plus_modu <- cluster_infomap(graphAD_plus, e.weights = graphAD_plus$mut_info_norm)

graphnoAD_plus_modu <- cluster_infomap(graphnoAD_plus, e.weights = graphnoAD_plus$mut_info_norm)

#Modules from plus networks

graphAD_plus_modules <- membership(graphAD_plus_modu)

graphnoAD_plus_modules <-  membership(graphnoAD_plus_modu)

# Assign the modules as attributes of the vertices

V(graphAD_plus)$modules <- graphAD_plus_modules

V(graphnoAD_plus)$modules <- graphnoAD_plus_modules

#Compare modularity between graphs, applying 
#variation of information "vi"
#normalized mutual information "nmi"
#split-join distance "split-join distance"
#Rand index "Rand index"
#adjusted Rand index "adjusted Rand index"

possible_algos <- c("vi", "nmi", "split.join", "rand", "adjusted.rand")

comparison_methods <- sapply(X = possible_algos, FUN = function(i){
  igraph::compare(comm1 = graphAD_plus_modules,
                  comm2 = graphnoAD_plus_modules,
                  method = i
  )
})

comparison_methods

#           vi           nmi    split.join          rand adjusted.rand 
#3.539291e+00  5.624657e-01  1.355000e+03  8.480256e-01  5.993397e-02 

# vi           nmi    split.join          rand adjusted.rand   <- 
#4.055584e+00  5.614082e-01  2.819000e+03  8.784207e-01  4.885023e-02

#END