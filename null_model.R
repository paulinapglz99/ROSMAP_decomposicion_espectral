#null_model.R
#This script builds two null models for each graph, by randomizing the networks while preserving the degree distribution. 
#Specifically, we rewired the edges of the original networks 1,000 times using a degree-preserving randomization algorithm. 
pacman::p_load("igraph", 
        "ggplot")

#Declare functions ----

# Function to calculate various metrics
calculate_metrics <- function(graph) {
  list(
    clustering = transitivity(graph, type = "global"),  # Global clustering coefficient
    avg_path_length = mean_distance(graph, directed = FALSE),  # Average path length
    modularity = modularity(cluster_infomap(graph)),  # Modularity of the network, using infomap
  )
}

#Function to randomize networks ----

randomize_graph <- function(graph, num_randomizations, niter = 1000) {
  random_graphs <- vector("list", num_randomizations)
  for (i in 1:num_randomizations) {
    random_graphs[[i]] <- rewire(graph, with = keeping_degseq(niter = niter))
  }
  return(random_graphs)
}

#Function to calculate metrics from ranzomized networks

calculate_metrics_random_graphs <- function(graph_list) {
  lapply(graph_list, calculate_metrics)
}

#Function to plot any metric -----

plot_metric_distribution <- function(random_metrics, observed_metric, metric_name, graph_type) {
  random_df <- data.frame(value = random_metrics, type = paste0("Random (", graph_type, ")")) # Df of random variables
  observed_df <- data.frame(value = observed_metric, type = paste0("Observed (", graph_type, ")"))   # Observed value
  combined_df <- rbind(random_df, observed_df)   # Combine data frames
  # Plot
  ggplot(random_df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.001, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "blue", lwd = 1) +
    geom_vline(aes(xintercept = observed_metric), color = "red", linetype = "dashed", size = 1.2) +
    labs(title = paste("Null Model for", metric_name, "(", graph_type, ")"), 
         x = metric_name, 
         y = "Density") + 
    theme_minimal()
}

#Read data ----

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

# Calculate metrics for your original networks
metrics_AD <- calculate_metrics(graphAD)
metrics_noAD <- calculate_metrics(graphnoAD)

# Number of random networks to generate
num_randomizations <- 1000

#Randomize

random_graphs_AD <- randomize_graph(graphAD, num_randomizations)
random_graphs_noAD <- randomize_graph(graphnoAD, num_randomizations)

# Calculate metrics for randomized graphs
metrics_AD_random <- calculate_metrics_random_graphs(random_graphs_AD)
metrics_noAD_random <- calculate_metrics_random_graphs(random_graphs_noAD)

# Extract metrics -----

# Extract clustering coefficient
clustering_AD_random <- sapply(metrics_AD_random, function(x) x$clustering)
clustering_noAD_random <- sapply(metrics_noAD_random, function(x) x$clustering)

# Plot clustering coefficient distribution
plot_metric_distribution(clustering_AD_random, metrics_AD$clustering, "Clustering Coefficient", "AD")
plot_metric_distribution(clustering_noAD_random, metrics_noAD$clustering, "Clustering Coefficient", "noAD")

# Extract average path length 

path_length_AD_random <- sapply(metrics_AD_random, function(x) x$avg_path_length)
path_length_noAD_random <- sapply(metrics_noAD_random, function(x) x$avg_path_length)

# Plot clustering average path length
plot_metric_distribution(path_length_AD_random, metrics_AD$avg_path_length, "Average path length", "AD")
plot_metric_distribution(path_length_noAD_random, metrics_noAD$avg_path_length, "Average path length", "noAD")

# Extract modularity for the random networks (AD and noAD)
modularity_AD_random <- sapply(metrics_AD_random, function(x) x$modularity)
modularity_noAD_random <- sapply(metrics_noAD_random, function(x) x$modularity)

# Plot clustering average path length
plot_metric_distribution(modularity_AD_random, metrics_AD$modularity, "Modularity", "AD")
plot_metric_distribution(modularity_noAD_random, metrics_noAD$modularity, "Modularity", "noAD")

# Confidence invtervals

#If calculate the 2.5 and 97.5 percentile for the randomized network metrics, we generate a 95% confidence interval. 
#If the observed metrics fall outside this interval, this will indicate that the observed metrics are not likely to be the product of chance.

# Intervalo de confianza al 95% para clustering coefficient
clustering_ci_AD <- quantile(clustering_AD_random, probs = c(0.025, 0.975))
clustering_ci_noAD <- quantile(clustering_noAD_random, probs = c(0.025, 0.975))

# Verificar si el valor observado cae fuera del intervalo de confianza
print(clustering_ci_AD)
print(metrics_AD$clustering)
print(clustering_ci_noAD)
print(metrics_noAD$clustering)

#END