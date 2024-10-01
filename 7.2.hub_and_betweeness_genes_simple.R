#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               'biomaRt', 
               "ggVennDiagram", 
               "clusterProfiler")

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

#Functions --- ---

# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define function to convert from ENSMBL to SYMBOL
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
  V(graph)$new_name <- sorted_graph_vnames_trad
  
  return(graph)
}

#Translate graphs

graph_AD <- translate_vertex_names(graphAD)

graphnoAD <- translate_vertex_names(graphnoAD)

##### HUB GENES ####

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Threshold
cutoff_degrees <- lapply(nodes_degree, function(deg) {
  quantile(deg, probs = 0.95)
})

#Hub nodes
hub_nodes <- mapply(function(graph, deg, cutoff) {
  # Extract nodes with degree >= 95th percentile
  V(graph)[deg >= cutoff]$name
}, graphs, nodes_degree, cutoff_degrees)

#What genes do both networks share?

shared_hub_genes <- intersect(hub_nodes[[1]], hub_nodes[[2]])

#What genes are in the AD network but not in the no AD network?

AD_hubs_notcontrol <- setdiff(hub_nodes[[1]], hub_nodes[[2]])
length(AD_hubs_notcontrol)

#Venn Diagram

ggVennDiagram(list(AD = hub_nodes[[1]], noAD =  hub_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

#Enrichment of these genes

AD_hubs_notcontrol_enrichment <- enrichGO(
  gene = AD_hubs_notcontrol,
  OrgDb = org.Hs.eg.db, 
  keyType = 'ENSEMBL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

#Plot enrichment of hubs

dotplot(AD_hubs_notcontrol_enrichment, showCategory=10)

cnetplot(AD_hubs_notcontrol_enrichment, circular = TRUE, colorEdge = TRUE) 

#Induced subgraph of hub genes

hub_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_hubs_notcontrol])

plot(hub_subgraph)

#Save graph 

write_graph(hub_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_hubs_induced_subgraph.graphml',
             format = "graphml")

##### HIGH BETWEENESS GENES ####

#Calculate degree of nodes
nodes_betwe <- sapply(X = graphs, FUN = betweenness)

#Threshold
cutoff_betwe <- lapply(nodes_betwe, function(x) {
  quantile(x, probs = 0.95)
})

#Hub nodes
highbe_nodes <- mapply(function(graph, x, cutoff) {
  # Extract nodes with degree >= 95th percentile
  V(graph)[x >= cutoff]$name
}, graphs, nodes_betwe, cutoff_betwe)

#What genes do both networks share?

shared_highbe_genes <- intersect(highbe_nodes[[1]], highbe_nodes[[2]])
shared_highbe_genes
#What genes are in the AD network but not in the no AD network?

AD_highbe_notcontrol <- setdiff(highbe_nodes[[1]], highbe_nodes[[2]])
length(AD_hubs_notcontrol)

#Venn Diagram

ggVennDiagram(list(AD = highbe_nodes[[1]], noAD =  highbe_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

#Enrichment of these genes

AD_highbe_notcontrol_enrichment <- enrichGO(
  gene = AD_highbe_notcontrol,
  OrgDb = org.Hs.eg.db, 
  keyType = 'ENSEMBL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

#Plot enrichment of hubs

dotplot(AD_highbe_notcontrol_enrichment, showCategory=10)

cnetplot(AD_highbe_notcontrol_enrichment, circular = TRUE, colorEdge = TRUE) 

#Induced subgraph of hub genes

highbeb_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_highbe_notcontrol])

plot(highbeb_subgraph)

#Save graph

write_graph(hub_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_highbe_induced_subgraph.graphml',
            format = "graphml")

#NEXT QUESTION IS What modules do these genes belong to? --- ---

#Save hub genes to explore them in the partitions

write(AD_hubs_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_hubs_notcontrolens.txt')

#Save high betweeness genes to explore them in the partitions

write(AD_highbe_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_ens.txt')

#END