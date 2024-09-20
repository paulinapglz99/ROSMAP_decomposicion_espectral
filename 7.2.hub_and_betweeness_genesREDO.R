#non sequitur
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
               "ClusterProfiler")

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

#Define functions --- --- 

#1. Function to convert gene names

# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

# Define function to change names of vertex for symbol names --- ---

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


#3. Create degree data frames

process_distribution <- function(degree_distribution) {
  # Convert to data frame and count frequencies
  distribution.df <- as.data.frame(table(degree_distribution$degree))
  # Rename the column
  distribution.df <- distribution.df %>% rename(degree = "Var1")
  # Order by degree
  distribution.df <- distribution.df[order(distribution.df$degree, decreasing = F), ]
  #Calculate mean and sd
  dis_mean <- mean(distribution.df$degree)
  dis_sd <- sd(distribution.df$degree)
  # Calculate the cumulative sum
  distribution.df$CumulativeDegree <- cumsum(distribution.df$Freq)
  # Logarithm of the cumulative sum
  distribution.df$logCumulativeDegree <- log10(distribution.df$CumulativeDegree)

  # Add threshold
  #simple_degree_threshold <- quantile(distribution.df$degree, probs = 0.95)
  logCumulativeDegree_threshold <- quantile(distribution.df$logCumulativeDegree, probs = 0.95)
  # Puedes devolver el umbral si es necesario, por ejemplo:
  return(list(degree_distribution = distribution.df, dis_mean= dis_mean,
              dis_sd = dis_sd, CumulativeDegree_threshold = CumulativeDegree_threshold))
}

########## HUB GENES ########## 

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Table of degree distribution

degree_distribution <- list()

for (i in 1:length(nodes_degree)) {
  degree_distribution[[i]] <- data.frame(gene = names(nodes_degree[[i]]), degree = nodes_degree[[i]])
}

#Calculate percentile 95 of genes with higher degree ---- ---

processed_distribution <- lapply(degree_distribution, FUN = process_distribution)

#Plot cummulative degree distribution --- ---

dis_AD <- processed_distribution[[1]]$degree_distribution
dis_AD$dx <- "AD"
dis_noAD <-  processed_distribution[[2]]$degree_distribution
dis_noAD$dx<- "no AD"

dis <- bind_rows(dis_AD, dis_noAD)
dis$degree <- as.numeric(dis$degree)

#Plot both distributions

degree_distr.p <- ggplot(dis, aes(x = degree, y = logCumulativeDegree, group = dx, color = dx)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = processed_distribution[[1]]$threshold, linetype = "dashed", color = "#291F1E", size = 0.5) +  # AD threshold
  geom_hline(yintercept = processed_distribution[[2]]$threshold, linetype = "dashed", color = "#291F1E", size = 0.5) +   # no AD threshold
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "logCumulativeDegree") +
  scale_x_continuous(breaks = seq(1, max(dis$degree))) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_light()

degree_distr.p    

#Save plot 

ggsave("d_distribution.png", plot = degree_distr.p, width = 15,
       height = 4, units = "in", dpi = 300)

#So, hub genes are

hub_AD <- dis_AD %>% filter(processed_distribution[[1]]$threshold)

####### HIGH BETWEENESS GENES #######



