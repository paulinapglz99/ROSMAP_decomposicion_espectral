# Example using igraph
library(igraph, 
        ggplot)

# Assume graphAD and graphnoAD are the adjacency matrices for the two networks
# Load or construct the networks

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

# Function to calculate various metrics
calculate_metrics <- function(graph) {
  list(
    clustering = transitivity(graph, type = "global"),  # Global clustering coefficient
    avg_path_length = mean_distance(graph, directed = FALSE),  # Average path length
    modularity = modularity(cluster_fast_greedy(graph)),  # Modularity of the network
    degree_distribution = degree_distribution(graph)  # Degree distribution
  )
}

# Calculate metrics for your original networks
metrics_AD <- calculate_metrics(graphAD)
metrics_noAD <- calculate_metrics(graphnoAD)

# Print metrics for reference
print(metrics_AD)
print(metrics_noAD)

#Create a Null Model (Randomize Networks)

# Number of random networks to generate
num_randomizations <- 1000

# Initialize lists to store metrics for random networks
metrics_AD_random <- vector("list", num_randomizations)
metrics_noAD_random <- vector("list", num_randomizations)

# Generate random networks and calculate metrics
for (i in 1:num_randomizations) {
  # Rewire the network while keeping the degree sequence intact
  random_graph_AD <- rewire(graphAD, with = keeping_degseq(niter = 1000))
  random_graph_noAD <- rewire(graphnoAD, with = keeping_degseq(niter = 1000))
  
  # Calculate metrics for the random networks
  metrics_AD_random[[i]] <- calculate_metrics(random_graph_AD)
  metrics_noAD_random[[i]] <- calculate_metrics(random_graph_noAD)
}

# Extract clustering coefficients for the random networks (AD and noAD)
clustering_AD_random <- sapply(metrics_AD_random, function(x) x$clustering)
clustering_noAD_random <- sapply(metrics_noAD_random, function(x) x$clustering)

# Extract average path lengths for the random networks (AD and noAD)
avg_path_length_AD_random <- sapply(metrics_AD_random, function(x) x$avg_path_length)
avg_path_length_noAD_random <- sapply(metrics_noAD_random, function(x) x$avg_path_length)

# Extract modularity for the random networks (AD and noAD)
modularity_AD_random <- sapply(metrics_AD_random, function(x) x$modularity)
modularity_noAD_random <- sapply(metrics_noAD_random, function(x) x$modularity)

# Create a data frame for graphAD random clustering coefficients
clustering_AD_df <- data.frame(clustering_coefficient = clustering_AD_random, 
                               type = "Random (AD)")

# Add the observed value to the data frame (for overlaying in the plot)
observed_AD_df <- data.frame(clustering_coefficient = metrics_AD$clustering, 
                             type = "Observed (AD)")

# Combine the random and observed values into one data frame
clustering_AD_combined_df <- rbind(clustering_AD_df, observed_AD_df)

# Create a similar data frame for graphnoAD
clustering_noAD_df <- data.frame(clustering_coefficient = clustering_noAD_random, 
                                 type = "Random (noAD)")

observed_noAD_df <- data.frame(clustering_coefficient = metrics_noAD$clustering, 
                               type = "Observed (noAD)")

clustering_noAD_combined_df <- rbind(clustering_noAD_df, observed_noAD_df)

# Plot distribution of clustering coefficients for AD graph

clustering_AD.p <- ggplot(clustering_AD_df, aes(x = clustering_coefficient)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "blue", lwd = 1) +  # Agrega una curva de densidad para mayor suavidad
  geom_vline(aes(xintercept = metrics_AD$clustering), color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Null Model for Clustering Coefficient (AD)", 
       x = "Clustering Coefficient", 
       y = "Density") + 
  theme_minimal()

# Plot for graphnoAD
clustering_noAD.p <- ggplot(clustering_noAD_df, aes(x = clustering_coefficient)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "blue", lwd = 1) +  # Agrega una curva de densidad para mayor suavidad
  geom_vline(aes(xintercept = metrics_AD$clustering), color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Null Model for Clustering Coefficient (AD)", 
       x = "Clustering Coefficient", 
       y = "Density") + 
  theme_minimal()


