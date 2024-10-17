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

##### Histograms #####

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

#Superponed

# Combinar ambos data.frames
combined_degree_freq <- rbind(ADdegree_freq, noADdegree_freq)

# Crear el gráfico combinado
log_log_combined <- ggplot(combined_degree_freq, aes(x = degree, y = Prob, color = dx)) +
  geom_point(size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  scale_x_log10() +  # Escala logarítmica en el eje x
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = " ", 
       subtitle = "a)") +
  geom_smooth(method = "lm", se = FALSE, linetype = "longdash", 
              size = 1.2, aes(color = dx), formula = y ~ x) +  # Ajuste lineal con colores diferenciados
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  ) +
  scale_color_manual(values = c("AD" = "blue", "No AD" = "red")) +
  theme(legend.position = "none" )

# Mostrar el gráfico
log_log_combined

ggsave(filename = "combined_degree_dis_loglog.jpg",
       plot = log_log_combined,
       width = 25,
       height = 25,
       units = "cm",
       dpi = 300,
)

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

#Cummulative degree

degree_distr_x.p <- ggplot(dis, aes(x = as.numeric(degree), y = as.numeric(CumulativeDegree), color = dx)) +
 #geom_point(size = 1.5) +
  geom_line(size = 2) +
  #geom_smooth(method = "loess", se = FALSE, size = 1.5) +  # Suavizado LOESS para suavizar la línea
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 15)) +  # Ajusta el eje X
  labs(subtitle = "b)", x = "Degree",
       y = "Cumulative Degree") +
  #scale_x_continuous(breaks = seq(1, max(dis$degree), by =4)) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Elimina las líneas de la cuadrícula principal
    #panel.grid.minor = element_blank()   # Elimina las líneas de la cuadrícula secundaria
  )

degree_distr_x.p


ggsave(filename = "accum_distr.jpg",
       plot = degree_distr_x.p,
       width = 50,
       height = 25,
       units = "cm",
       dpi = 300,
)


#Final panel

panel <- grid.arrange(log_log_combined, degree_distr_x.p, ncol = 2)

ggsave(filename = "panel_loglog_accum_distr.jpg",
       plot = panel,
       width = 50,
       height = 25,
       units = "cm",
       dpi = 300,
)

#END