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

###### log-log graph ######

log_log_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
 scale_x_log10() +  # Escala logarítmica en el eje x
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", 
               size = 1.2, formula = y ~ x) +  # Ajuste lineal con estilo modificado
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )


log_log_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
 scale_x_log10() + 
 scale_y_log10() + 
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "No AD network \ndegree distribution") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", 
               size = 1.2, formula = y ~ x) +  # Ajuste lineal con estilo modificado
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

log_grid <- grid.arrange(log_log_AD, log_log_noAD, ncol =2)

# ggsave(filename = "degree_dis_loglog.jpg",
#        plot = log_grid,
#        width = 25,
#        height = 25,
#        units = "cm",
#        dpi = 300,
# )


###########  k vs P(k) distribution ########

k_pk_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(p(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(p(k)), 
       title = "", 
       subtitle = "No AD network \ndegree distribution") +
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_grid <- grid.arrange(k_pk_AD, k_pk_noAD, ncol =2)

###### k vs logp(k) #####

k_logpk_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_logpk_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "No AD network \ndegree distribution") +
  scale_y_log10() +  # Escala logarítmica en el eje y
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_grid <- grid.arrange(k_logpk_AD, k_logpk_noAD, ncol =2)

#END