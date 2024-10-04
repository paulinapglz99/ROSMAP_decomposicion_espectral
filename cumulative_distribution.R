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

#Calculate degree of nodes
nodes_degree <- sapply(X = graphLists, FUN = degree)

#Table of degree distribution
degree_distribution <- list()

for (i in 1:length(nodes_degree)) {
  degree_distribution[[i]] <- data.frame(gene = names(nodes_degree[[i]]), degree = nodes_degree[[i]])
}

#
processed_distribution <- lapply(degree_distribution, FUN = process_distribution)

dis_AD <- processed_distribution[[1]]$degree_distribution
dis_AD$dx <- "AD"
dis_AD$degree <- as.numeric(dis_AD$degree)
max(dis_AD$degree)
#[1] 224

dis_noAD <-  processed_distribution[[2]]$degree_distribution
dis_noAD$dx<- "no AD"
dis_noAD$degree <- as.numeric(dis_noAD$degree)
max(dis_noAD$degree)

#[1] 223

dis<-  bind_rows(dis_AD, dis_noAD)
max(as.numeric(dis$degree))


degree_distr.p <- ggplot(dis, aes(x = degree, y = logCumulativeDegree, color = dx)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 5)) +  # Ajusta el eje X
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "logCumulativeDegree") +
  #scale_x_continuous(breaks = seq(1, max(dis$degree), by =4)) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_minimal()

degree_distr.p

#Save plot

ggsave(filename = "bothdx_cumulative_degree_distributions_coexpression_NIAReagan_histogram.png",
       plot = degree_distr.p,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300,
)

#

degree_distr_x.p <- ggplot(dis, aes(x = degree, y = CumulativeDegree, color = dx)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 5)) +  # Ajusta el eje X
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "Cumulative Degree") +
  #scale_x_continuous(breaks = seq(1, max(dis$degree), by =4)) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_minimal()

degree_distr_x.p


ggsave(filename = "bothdx_nolog_cumulative_degree_distributions_coexpression_NIAReagan_histogram.png",
       plot = degree_distr_x.p,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300)

#END