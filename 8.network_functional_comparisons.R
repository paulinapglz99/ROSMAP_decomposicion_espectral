#
#8.network_functional_comparisons.R

#This script makes the functional modular comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#Modules will be compared in terms of the associated biological functions identified by an enrichment analysis.
#This comparison involves answering three complementary questions: 

#1.How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?
#2.How similar are the modules found in each network, in terms of the sets of associated biological functions?
#3.In how many modules is represented each biological process?

#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse",
               "clusterProfiler", 
               "gridExtra", 
               "biomaRt")

library("org.Hs.eg.db", character.only = TRUE)

#Set seed for modularity algorithm --- ---

set.seed(10)

# Define functions --- ---

# Connecting to the Ensembl database through biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
}

# Define function for translating vertex names 

translate_vertex_names <- function(graph) {
  graph_vnames <- V(graph)$name   # Extract vertex names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)   # Translate names
  # Replace the missing values in the column 'external_gene_name' with the values of 'ensembl_gene_id'.
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  # Sort graph_vnames_trad according to the order of graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
  # Assign the new names to the network vertices.
  V(graph)$name <- sorted_graph_vnames_trad
  return(graph)
}

#Define similarity of Enriched Processes, Jaccard Index function --- ---

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

# Crear una función para calcular el Jaccard Index
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

#Get data --- --- 

graphAD <- read_graph(file =  '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
                        format = 'graphml')

universe <- scan(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/universe.txt", what = character())

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Comparison of biological function sets associated to the overall network --- ---

#Enrichment of full graphs

#For AD graph

enrichment_fullnet_AD <- enrichGO(gene = V(graphAD)$name,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = 'ENSEMBL',
                                  readable = TRUE,
                                  universe = universe, 
                                  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

#Gene concept network enrichment

enrichment_fullnet_AD_cnet <- cnetplot(enrichment_fullnet_AD, showCategory= 5)

#Dotplot enrichment

enrichment_fullnet_AD_dot <- dotplot(enrichment_fullnet_AD)

#Save plots

ggsave("enrichment_fullnet_AD_cnet.png", enrichment_fullnet_AD_cnet, width = 15, height = 8)

#For noAD graph

enrichment_fullnet_noAD <- enrichGO(gene = V(graphnoAD)$name,
                                    OrgDb = "org.Hs.eg.db", 
                                    keyType = 'ENSEMBL',
                                    readable = TRUE,
                                    universe = universe, 
                                    ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

#Gene concept network enrichment

enrichment_fullnet_noAD_cnet <- cnetplot(enrichment_fullnet_noAD, showCategory= 5)

#Dotplot enrichment

enrichment_fullnet_noAD_dot <- dotplot(enrichment_fullnet_noAD)

#Comparison of enrichments --- ---

OverallEnrichedProcessJ <- jaccard_simplex(names(enrichment_fullnet_AD@geneSets), names(enrichment_fullnet_noAD@geneSets))
#[1] 0.8042686

#This answers the question 2. How similar are the modules found in each network, in terms of the sets of associated biological functions?
  
#Similarity of modules in terms of associated biological functions --- ---

#Modularity algorithm 

modularities <- sapply(X = graphLists, FUN = cluster_infomap)

#Calculate Q score 

Qscore <- sapply(modularities, FUN = modularity)

#Calculate clustering coefficient

clus_coe <- sapply(graphLists, FUN = transitivity)

#Split lists of nodes by module

nodes_membership <- sapply(modularities, FUN = membership)

#Create df
nodes_membership_AD.df <- data.frame(ensembl_gene_id = names(nodes_membership$graphAD),  membership = nodes_membership$graphAD)

nodes_membership_noAD.df <- data.frame(ensembl_gene_id = names(nodes_membership$graphnoAD),  membership = nodes_membership$graphnoAD)

# Assign membership to the nodes of each network
for (i in seq_along(graphLists)) {
  V(graphLists[[i]])$community <- nodes_membership[[i]]
}

#Save graphs  --- ----

# Loop through each graph and export it with community membership as GraphML
for (i in seq_along(graphLists)) {
  # Generate filename for each graph (e.g., graphAD.graphml, graphnoAD.graphml)
  graph_name <- ifelse(i == 1, "graphAD", "graphnoAD")
  filename <- paste0(graph_name, "nodes_membership.graphml")
  
  # Export graph to GraphML format
  write_graph(graphLists[[i]], file = filename, format = "graphml")
}

# Extract list of nodes by community for each graph

nodes_by_community_list <- lapply(seq_along(nodes_membership), function(i) {
  split(V(graphLists[[i]])$name, nodes_membership[[i]])
})

# Extraer los genes de cada módulo
modules_AD <- split(nodes_membership_AD.df$ensembl_gene_id, nodes_membership_AD.df$membership)
modules_noAD <- split(nodes_membership_noAD.df$ensembl_gene_id, nodes_membership_noAD.df$membership)

# Crear la matriz de similitud
similarity_matrix <- matrix(0, nrow = length(modules_AD), ncol = length(modules_noAD))

# Rellenar la matriz con los índices de similitud Jaccard
for (i in seq_along(modules_AD)) {
  for (j in seq_along(modules_noAD)) {
    similarity_matrix[i, j] <- jaccard_index(modules_AD[[i]], modules_noAD[[j]])
  }
}

# Nombrar filas y columnas según los módulos
rownames(similarity_matrix) <- paste0("Module_AD_", seq_along(modules_AD))
colnames(similarity_matrix) <- paste0("Module_noAD_", seq_along(modules_noAD))

# Visualizar la matriz de similitud
similarity_matrix
dim(similarity_matrix)
#[1] 65 71

#Is there any modules with perfect similitude?

length(similarity_matrix[similarity_matrix == 1])
#[1] 10

# Buscar posiciones de valores iguales a 1
pairs_with_one <- which(similarity_matrix == 1, arr.ind = TRUE)

# Crear dataframe con los nombres de los módulos correspondientes
module_pairs <- data.frame(
  module_AD = rownames(similarity_matrix)[pairs_with_one[, "row"]],
  module_noAD = colnames(similarity_matrix)[pairs_with_one[, "col"]]
)

# Imprimir el dataframe corregido
print(module_pairs)

#In which modules are found our hub genes? --- ---

#Get hub genes 

hub_genes <- scan("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_hubs_notcontrolens.txt",  what = character())

hub_genes.x <- nodes_membership_AD.df %>% filter(ensembl_gene_id %in% hub_genes)

hub_genes.xy <- convert_ens_to_symbol(hub_genes.x)

hub_genes.x <- hub_genes.x %>% left_join(hub_genes.xy, by ="ensembl_gene_id")

#In which modules are found our high betweeness genes?--- ---

high_be_genes <-  scan(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/high_be_genes_en_AD_no_noAD_ens.txt",   what = character())

high_be_genes.x <- nodes_membership_AD.df %>% filter(ensembl_gene_id %in% high_be_genes)

high_be_genes.xy <- convert_ens_to_symbol(high_be_genes.x)

high_be_genes.x <- high_be_genes.x %>% left_join(high_be_genes.xy, by ="ensembl_gene_id")

#Enrichment to all modules of a network --- ----

#I'd like to add names to lists of graphs again

names(nodes_by_community_list)[1] <- "graphAD"
names(nodes_by_community_list)[2] <- "graphnoAD"

#Extract list of nodes by community 

nodes_by_community_AD <- nodes_by_community_list[["graphAD"]]

nodes_by_community_noAD <- nodes_by_community_list[["graphnoAD"]]

#Define function that performs enrichment

enricher <- function(nodes_by_community) {
  enrichGO_result <- enrichGO(gene = nodes_by_community,
                              OrgDb = org.Hs.eg.db, 
                              universe = universe , 
                              keyType = 'ENSEMBL',
                              readable = TRUE,
                              ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)
  return(enrichGO_result)
}

#These functions may print a " --> No gene can be mapped...." on the console. 
#This happens when there are lists of genes that cannot be enriched. It will still generate a result

# Define a function to replace NULL by an empty S4 object "enrichResult" to handle NULLs in the enrichment lists

replace_null <- function(x) {
  if (is.null(x)) {
    return(new("enrichResult"))
  } else {
    return(x)
  }
}

#Enrichment of graphAD modules

enriched_results_AD <- lapply(nodes_by_community_AD, enricher) #slow
enriched_results_AD <- lapply(enriched_results_AD, replace_null) # Exchange NULLs for empty S4 objects

#Enrichment of graphAD modules
enriched_results_noAD <- lapply(nodes_by_community_noAD, enricher) #slow
enriched_results_noAD <- lapply(enriched_results_noAD, replace_null) # Exchange NULLs for empty S4 objects

#Functional module comparison between two graphs --- ---

#1.How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?

#Similarity of modules in terms of associated biological functions

#It is possible to compare each module of the network of interest to the modules
#of another network, in terms of their associated biological function

# Create an empty array to store the results of the Jaccard index.
num_modules_AD <- length(enriched_results_AD)
num_modules_noAD <- length(enriched_results_noAD)

similarity_matrix <- matrix(NA, nrow = num_modules_AD, ncol = num_modules_noAD)

# Assign row and column names
rownames(similarity_matrix) <- paste("AD", 1:num_modules_AD, sep = "_")
colnames(similarity_matrix) <- paste("noAD", 1:num_modules_noAD, sep = "_")

#Iterate to obtain similarity matrix 

for (i in 1:num_modules_AD) {
  for (j in 1:num_modules_noAD) {
    # Acceder a los geneSets de las redes AD y no AD
    geneSets_AD <- names(enriched_results_AD[[i]]@geneSets)
    geneSets_noAD <- names(enriched_results_noAD[[j]]@geneSets)
    similarity_matrix[i, j] <- jaccard_simplex(geneSets_AD, geneSets_noAD)
  }
}

head(similarity_matrix)

#Find the highest similarity of the matrix  
max(similarity_matrix, na.rm = TRUE)

#Network with this matrix

sim_node.g <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = TRUE)

#Number of modules in networks that have Jaccard index J=1 with a module of the Main network --- ---

# Function to count the number of columns with value equal to 1
count_ones <- function(row) {
  sum(row == 1)
}

# Apply the function to each row of the matrix
ones_count <- apply(similarity_matrix, 1, count_ones)

# Crear una tabla con los resultados
num_equal_nodes <- data.frame(module_number = rownames(similarity_matrix), modules_with_jaccardindex1 = ones_count)

num_equal_nodes.x <- num_equal_nodes %>% filter(modules_with_jaccardindex1 ==1)
dim(num_equal_nodes.x)
#[1] 15  2

#This answers the question of 1. How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?
  
#Number of modules associated to a given biological function --- ---

#3.In how many modules is represented each biological process?

# Define a function to extract the biological processes (GO terms) from an enrichment result
extract_bp_terms <- function(enrich_result) {
  if (is.null(enrich_result@result)) return(character(0)) # Handle empty enrichment results
  bp_terms <- enrich_result@result$Description # Extract Biological Process (BP) terms
  return(bp_terms)
}

# Combine both AD and noAD enrichment results
combined_enrichment_results <- c(enriched_results_AD, enriched_results_noAD)

# Create a list to store the modules for each biological process
biological_process_modules <- list()

# Loop through each module's enrichment result
for (i in seq_along(combined_enrichment_results)) {
  # Extract the biological process terms for the current module
  bp_terms <- extract_bp_terms(combined_enrichment_results[[i]])
  
  # Count the occurrence of each biological process and track the module it appears in
  for (bp in bp_terms) {
    if (!is.null(biological_process_modules[[bp]])) {
      biological_process_modules[[bp]]$modules <- c(biological_process_modules[[bp]]$modules, i)
      biological_process_modules[[bp]]$count <- biological_process_modules[[bp]]$count + 1
    } else {
      biological_process_modules[[bp]] <- list(modules = i, count = 1)
    }
  }
}

# Convert the result to a data frame for easier viewing
biological_process_summary <- data.frame(
  Biological_Process = names(biological_process_modules),
  Module_Count = sapply(biological_process_modules, function(x) x$count),
  Modules = sapply(biological_process_modules, function(x) paste(x$modules, collapse = ", "))
)

# Display the result
print(biological_process_summary)



