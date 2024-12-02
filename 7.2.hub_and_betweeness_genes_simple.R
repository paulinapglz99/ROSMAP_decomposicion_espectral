#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               "ggVennDiagram", 
               "clusterProfiler", 
               "ggraph", 
               "biomaRt")

library("org.Hs.eg.db", character.only = TRUE)

#Function to translate genes --- ---
# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

translate_vertex_names <- function(graph) {
  #Extract vertex names
  graph_vnames <- V(graph)$name
  
  #Translate names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)
  
  #Replace the missing values in the column 'external_gene_name' with the values of 'ensembl_gene_id'.
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  
  # Sort graph_vnames_trad according to the order of graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
  #Assign the new names to the network vertices.
  V(graph)$name_trad <- sorted_graph_vnames_trad
  
  return(graph)
}

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

##### HUB GENES ####

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Threshold
cutoff_degrees <- lapply(nodes_degree, function(deg) {
  quantile(deg, probs = 0.90)
})

#Hub nodes

hub_nodes <- list()
for (name in names(graphs)) {
  deg <- nodes_degree[[name]]
  cutoff <- cutoff_degrees[[name]]
  hub_nodes[[name]] <- V(graphs[[name]])[deg >= cutoff]$name
}

#What genes do both networks share?

shared_hub_genes <- intersect(hub_nodes[[1]], hub_nodes[[2]])

#What genes are in the AD network but not in the no AD network?

AD_hubs_notcontrol <- setdiff(hub_nodes[[1]], hub_nodes[[2]])
length(AD_hubs_notcontrol)

AD_hubs_notcontrol.df1 <- data.frame(gene = names(nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                        AD_hubs_notcontrol]),
                                    degree = nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                    AD_hubs_notcontrol]) 

AD_hubs_notcontrol.df <- convert_ens_to_symbol(AD_hubs_notcontrol.df1$gene)


AD_hubs_notcontrol.df <- AD_hubs_notcontrol.df %>%
  left_join(AD_hubs_notcontrol.df1, by = c("ensembl_gene_id" = "gene")) %>% arrange(desc(degree))

#Venn Diagram

hubsvenn <- ggVennDiagram(list(AD = hub_nodes[[1]], C =  hub_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "a) Shared Hub genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
hubsvenn

#Enrichment of these genes

AD_hubs_notcontrol_enrichment <- enrichGO(
  gene = AD_hubs_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  universe = universe,
  ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

AD_hubs_notcontrol_enrichment.df <- as.data.frame(AD_hubs_notcontrol_enrichment)

#Plot enrichment of hubs

cnetplot(AD_hubs_notcontrol_enrichment, circular = F, colorEdge = TRUE) 

#Induced subgraph of hub genes AD_hubs_notcontrol

V(graphAD)$degree <- nodes_degree[["graphAD"]]

hub_subgraph_AD <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% hub_nodes[["graphAD"]]])

plot(hub_subgraph_AD)

hub_subgraph_AD.p <- ggraph(hub_subgraph_AD, layout = 'kk') + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = degree), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Degree",low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("a)") +
  theme_void() 
hub_subgraph_AD.p 

#Induced subgraph of hub genes AD_hubs_notcontrol

hub_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_hubs_notcontrol])

hub_subgraph.p <- ggraph(hub_subgraph, layout = 'kk') + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = degree), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Degree", low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("c)") +
  theme_void() 
hub_subgraph.p

#Save graph 

write_graph(hub_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_hubs_induced_subgraph.graphml',
             format = "graphml")

#Enrichment with KEGG

kk <- enrichKEGG(gene = AD_hubs_notcontrol,
                 organism = 'hsa',
                 keyType = 'ENSEMBL',
                 pvalueCutoff = 0.05, 
                 universe = universe)

##### HIGH BETWEENESS GENES ####

#Calculate degree of nodes
nodes_betwe <- sapply(X = graphs, FUN = betweenness)

#Threshold
cutoff_betwe <- lapply(nodes_betwe, function(x) {
  quantile(x, probs = 0.88)
})

highbe_nodes <- list()
for (name in names(graphs)) {
  betwe <- nodes_betwe[[name]]
  cutoff <- cutoff_betwe[[name]]
  highbe_nodes[[name]] <- V(graphs[[name]])[betwe >= cutoff]$name
}

#What genes do both networks share?

shared_highbe_genes <- intersect(highbe_nodes[[1]], highbe_nodes[[2]])
shared_highbe_genes
#What genes are in the AD network but not in the no AD network?

AD_highbe_notcontrol <- setdiff(highbe_nodes[[1]], highbe_nodes[[2]])
length(AD_highbe_notcontrol)

#Venn Diagram

highbevenn <- ggVennDiagram(list(AD = highbe_nodes[[1]], C =  highbe_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "b) Shared high betweeness genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

#Enrichment of these genes

AD_highbe_notcontrol_enrichment <- enrichGO(
  gene = AD_highbe_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

AD_highbe_notcontrol_enrichment.df <- as.data.frame(AD_highbe_notcontrol_enrichment)

vroom::vroom_write(AD_highbe_notcontrol_enrichment.df,
                   "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_enrichment.csv")

#Plot enrichment of hubs

AD_highbe_notcontrol_enrichment.p <- cnetplot(AD_highbe_notcontrol_enrichment, circular = TRUE, colorEdge = TRUE) +
  guides(edge_color = "none") +
  ggtitle("e)") 

ggsave("AD_highbe_notcontrol_enrichment.jpg",
       plot = AD_highbe_notcontrol_enrichment.p, 
       device = "jpg",
       width = 15, #55
       height = 10, #30
       units = "in", 
       dpi = 300
)

#Induced subgraph of hub genes

V(graphAD)$betweenness <- nodes_betwe[["graphAD"]]

highbe_subgraph_AD <- induced_subgraph(graph = graphAD, 
                                    vids = V(graphAD)[name %in% highbe_nodes[["graphAD"]]])

highbe_subgraph_AD.p <- ggraph(highbe_subgraph_AD, layout = 'kk', maxiter = 1000) + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = betweenness), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Betweenness centrality", low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("c)") +
  theme_void() 

highbe_subgraph_AD.p
#Genes only in AD graph

highbeb_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_highbe_notcontrol])

plot(highbeb_subgraph)

#Save graph

write_graph(highbeb_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_highbe_induced_subgraph.graphml',
            format = "graphml")

#Only in AD

highbeb_subgraph.p <- ggraph(highbeb_subgraph, layout = 'kk') + 
  geom_edge_link(aes(alpha = 0.8), color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = betweenness), size = 5) +  # Dibujar nodos sin leyenda
  scale_color_gradient("Betweenness", low = "lightblue", high = "red") +  # Eliminar la leyenda del color
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("d)") +
  theme_void() +
  theme(legend.title = element_text(hjust = 0.9))

#Plot graphs and Venn diagrams in a grid

lay <- rbind(c(3,1,1,2,2),  # Primera fila con tres gráficos
             c(5,4,4,2,2))  # Segunda fila con dos gráficos más, debajo del gráfico central

 grid <- grid.arrange(hub_subgraph.p,
                      highbeb_subgraph.p, 
                      hubsvenn,
                      AD_highbe_notcontrol_enrichment.p, 
                      highbevenn,
                      layout_matrix = lay)

 ggsave("grid_hubs_and_highbe1.jpg",
        plot = grid, 
        device = "jpg",
        width = 18, #55
        height = 12, #30
        units = "in", 
        dpi = 300)
 
 #NEXT QUESTION IS What modules do these genes belong to? --- ---

#Save hub genes to explore them in the partitions

write(AD_hubs_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_hubs_notcontrolens.txt')

#Save high betweeness genes to explore them in the partitions

write(AD_highbe_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_ens.txt')

#END