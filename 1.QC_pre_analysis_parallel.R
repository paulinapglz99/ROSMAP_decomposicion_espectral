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

#[1] "frontal cortex"                 "temporal cortex"  "dorsolateral prefrontal cortex" "Head of caudate nucleus"       
#[5] "posterior cingulate cortex"    

counts_FC <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FC/ROSMAP_RNAseq_rawcounts_FC.rds") %>% as.data.frame()
counts_TC <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/TC/ROSMAP_RNAseq_rawcounts_TC.rds") %>% as.data.frame()
counts_DLFPC <-readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds") %>% as.data.frame()
counts_HCN <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/HCN/ROSMAP_RNAseq_rawcounts_HCN.rds") %>% as.data.frame()
counts_PCC <- readRDS(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.rds") %>% as.data.frame()

counts <- list(counts_FC = counts_FC, 
               counts_TC = counts_TC,
               counts_DLFPC = counts_DLFPC,
               counts_HCN = counts_HCN, 
               counts_PCC = counts_PCC
               )

lapply(counts, dim)

#Obtain factors from metadata

metadata_FC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FC/RNA_seq_metadata_FC.txt")
metadata_TC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/TC/RNA_seq_metadata_TC.txt")
metadata_DLFPC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
metadata_HCN <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN/RNA_seq_metadata_HCN.txt")
metadata_PCC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.txt")

metadata <- list(metadata_FC = metadata_FC, 
                 metadata_TC = metadata_TC,
                 metadata_DLFPC = metadata_DLFPC, 
                 metadata_HCN = metadata_HCN, 
                 metadata_PCC = metadata_PCC)

lapply(metadata, dim)

#Change name of the first column

counts <- lapply(counts, function(counts) {
  names(counts)[1] <- "feature"
  return(counts)
})

#Trim STAR alignment stats

counts <- lapply(counts, function(x) x[-c(1:4), ])
lapply(counts, dim)

# Delete the "_PAR_Y" elements

counts <- lapply(counts, function(x) x[-grep("_PAR_Y", x$feature), ])
lapply(counts, dim)

#Prepare feature column
counts <- lapply(counts, function(x) {
  x <- x %>% mutate(feature = str_remove(feature, "\\..*$"))
  return(x)
})

#Create sub metadata

metadata <- lapply(metadata, function(x) {
  x %>% filter(!is.na(sequencingBatch))
})

lapply(metadata, dim)

#Filter counts to have only specimenes with metadata, just to be sure

counts <- mapply(function(counts_matrix, metadata_df) {
  specimenIDs <- metadata_df$specimenID
  counts_matrix %>%
    dplyr::select(c(1, any_of(specimenIDs)))
}, counts, metadata, SIMPLIFY = FALSE)

lapply(counts, dim)

# Delete duplicates

del_dupl <- function(counts) {
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
counts <- lapply(counts, del_dupl)

# PCA --- ---

#Transpose matrices 

#Function to build matrices
pca_matrix.f <- function(x) {
    x %>% 
    column_to_rownames("feature") %>% 
    as.matrix() %>% 
    t() # transpose the matrix so that rows = samples and columns = variables, this because dots in the PCA scatterplot will be the ones in the rows
}

#Build matrices

pca_matrix.l <- lapply(counts, pca_matrix.f) #rows must be = samples and columns = variables

# Look at the first 10 rows and first 5 columns of the matrix
lapply(pca_matrix.l, function(x) x[1:10, 1:10])

# PCA

pca.l <- lapply(pca_matrix.l, function(x) prcomp(x, retx = TRUE, center = TRUE, scale. = FALSE)) #slow

#PCA to table
pca_df.l <- lapply(pca, function(x) x$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID'))

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

var_table.df <- lapply(pca_list, function(pca) {
  # Calculate the variance of each component
  variance_percentage <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  # Calculate the cumulative variance
  cumulative_percentage <- cumsum(variance_percentage)
  # Create a dataframe with the information
  data.frame(
    PC = 1:length(pca$sdev),
    Variance_Percentage = variance_percentage,
    Cumulative_Percentage = cumulative_percentage
  )
})

#Create elbow plots

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

#PCA scatterplot coloured by library batch function

scatterplot_librarybatch.f <- function(pca_df, metadata, variance_table) {
  pca_df %>%
    ggplot() +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$sequencingBatch)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = "PCA Scatterplot coloured by library batch",
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
    theme_minimal()
}

#PCA scatterplot coloured by NIA Reagan function

scatterplot_NIAReagan.f <- function(pca_df, metadata_df, var_table.df) {
  pca_df %>%
    ggplot() +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$dicho_NIA_Reagan)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = "PCA Scatterplot coloured by library batch",
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table.df$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table.df$Variance_Percentage[2]), "%)")) +
    theme_minimal()
}

#Scatterplot coloured by sequencing batch

pca_seqbatch <- mapply(function(pca_df, metadata_df) {
  scatterplot_librarybatch.f(pca_df, metadata, var_table.df)
}, pca_df.l, metadata, SIMPLIFY = FALSE)


#Scatterplot coloured by dx

pca_NIAReagan <- mapply(function(pca_df, metadata_df) {
  scatterplot_NIAReagan.f(pca_df, metadata, var_table.df)
}, pca_df.l, metadata, SIMPLIFY = FALSE)
