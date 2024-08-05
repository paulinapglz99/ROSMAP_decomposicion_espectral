#
#1.QC_pre_analysis.R
#Script to inspect and pre-process RNA-seq data counts 
#paulinapglz.99@gmail.com

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

###Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2", 
               "ggfortify", 
               "gridExtra", 
               "edgeR")

#Get data --- --- 

counts_DLFPC_ROSMAP <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts_DLPFC_MSBB <-
counts_DLFPC_Mayo <- 

counts <- list(counts_DLFPC = counts_DLFPC, 
               counts_FC = counts_FC, 
               counts_HCN = counts_HCN, 
               counts_PCC = counts_PCC, 
               counts_TC = counts_TC)

lapply(counts, dim)

#Obtain factors from metadata

metadata_DLFPC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt") %>% as.data.frame()
metadata_FC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FC/RNA_seq_metadata_FC.txt") %>% as.data.frame()
metadata_HCN <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN/RNA_seq_metadata_HCN.txt") %>% as.data.frame()
metadata_PCC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.txt") %>% as.data.frame()

metadata <- list(metadata_DLFPC = metadata_DLFPC, 
                 metadata_FC = metadata_FC, 
                 metadata_HCN = metadata_HCN, 
                 metadata_PCC = metadata_PCC)

lapply(metadata, dim)

#Trim STAR alignment stats

counts <- lapply(counts, function(x) x[-c(1:4), ])
lapply(counts, dim)

# Delete the "_PAR_Y" elements

counts <- lapply(counts, function(x) -grep("_PAR_Y", x$feature))
lapply(counts, dim)

#Prepare feature column

counts <- lapply(counts, function(x) {
  x %>% mutate(feature = str_remove(feature, "\\..*$"), .before = 1)
})

#Create sub metadata

metadata <- lapply(metadata, function(x) {
  x %>% filter(!is.na(sequencingBatch))
})

#Filter counts to have only specimenes with metadata

counts <- lapply(counts, function(x) {
  x %>%
    dplyr::select(c(1, any_of(metadata$specimenID)))
})

# Delete duplicates

process_counts <- function(counts) {
  # Verificar genes repetidos
  repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)
  
  # Ver las filas duplicadas
  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  
  # Ordenar y calcular la mediana de los valores duplicados
  repeated_rows <- repeated_rows[order(repeated_rows$feature),]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), median, na.rm = TRUE))
  
  # Eliminar las filas duplicadas y agregar las filas con la mediana calculada
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  
  # Convertir las columnas seleccionadas a enteros
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  # Establecer los nombres de las filas
  rownames(counts) <- counts$feature
  
  return(counts)
}

# Aplicar la función a cada dataframe en la lista
counts <- lapply(counts, process_counts)

# PCA --- ---

#Transpose matrices 

matrices <- 

#Build matrix
  
pca_matrix.f <- function(x) {
    x %>% 
    column_to_rownames("feature") %>% 
    as.matrix() %>% 
    t()
}

pca_matrix.l <- lapply(df_list, convert_to_pca_matrix)
# Look at the first 10 rows and first 5 columns of the matrix
lapply(pca_matrix.l, function(x) x[1:10, 1:10])

# PCA

pca <- lapply(pca_matrix.l, function(x) prcomp(x, retx = TRUE, center = TRUE, scale. = FALSE)) #slow

#PCA to table
pca_df <- lapply(pca, function(x) x$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID'))

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

var_table.f <- function(pca) {
  # Calcular la varianza de cada componente
  variance_percentage <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  # Calcular la varianza acumulada
  cumulative_percentage <- cumsum(variance_percentage)
  # Crear un dataframe con la información
  data.frame(
    PC = 1:length(pca$sdev),
    Variance_Percentage = variance_percentage,
    Cumulative_Percentage = cumulative_percentage
  )
}

# Aplicar la función a cada objeto PCA en la lista
var_table.df <- lapply(pca_list, var_table.f)

elbow_plots.l <- lapply(var_table.df, function(x){
  ggplot(x, aes(x = PC, y = Variance_Percentage)) +
    geom_bar(stat = 'identity', fill = x$PC, color = 'black', position = 'dodge') +
    labs(title = 'Scree plot',
         subtitle = 'filtered data',
         x = 'Principal Components',
         y = 'Variance percentage') +
    scale_x_discrete() +
    theme_minimal()
})

#

scatterplot.l <- 

PC1_PC2_librarybatch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(metadata$sequencingBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = paste("PC1 vs PC2"),
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()
PC1_PC2_librarybatch