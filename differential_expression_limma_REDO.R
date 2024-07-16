#differential_expression_limma.R
#This script performs differential expression with limma-voom 

#Libraries --- ---

# BiocManager::install("limma")
# BiocManager::install("Glimma")
# BiocManager::install("edgeR")
# BiocManager::install("EnhancedVolcano")

pacman::p_load("tidyverse",
               "limma",
               "Glimma", 
               "edgeR", 
               "stringr", 
               "pheatmap",
               "edgeR",
               "EnhancedVolcano")

#Get data --- ---

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
counts$feature <- str_remove(counts$feature, "\\..*$")

# Combine and summarize duplicate features, then replace the original rows
# counts <- counts %>%
#   group_by(feature) %>% #group the data by the feature column.
#   summarize(across(everything(), #This combines the duplicate rows by calculating the median of their values.
#                    median, na.rm = TRUE),
#             .groups = 'drop')    #This ensures that the result does not maintain the grouping.
# dim(counts)
# #[1] 60558  1142

# Verify repeated genes in 'feature' column
repeated_values <- counts %>%
  group_by(feature) %>%
  filter(n() > 1) %>%
  distinct(feature) %>%
  pull(feature)
length(repeated_values)

# Ver las filas duplicadas
repeated_rows <- counts[counts$feature %in% repeated_values, ]
dim(repeated_rows)

repeated_rows <- repeated_rows[order(repeated_rows$feature),]

repeated_rows <- repeated_rows %>%
  group_by(feature) %>%
  summarize(across(everything(), median, na.rm = TRUE))
dim(repeated_rows)

#Delete rows in the matrix

counts <- counts %>% filter(!feature  %in% repeated_rows$feature)
dim(counts)
counts <- bind_rows(counts, repeated_rows)
dim(counts)
#[1] 60558  1142

rownames(counts) <- counts$feature
#  Convert selected columns to integers
counts <- counts %>% mutate(across(-feature, as.integer))

# Get metadata
metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   42

table(metadata$dicho_NIA_reagan, useNA = "ifany")
#  0    1 <NA> 
#  307  573  261  

#Filter to have only Samples with dicho_NIA_Reagan metadata

metadata <- metadata %>% filter(!is.na(dicho_NIA_reagan))

counts <- counts %>% dplyr::select(one_of(metadata$specimenID))
dim(counts)
#[1] 60558   793

#QC --- ---

# Library sizes

#Library sizes are a technical bias, to be corrected.
#This plot should look not that variant
barplot(dge$samples$lib.size)

# Create dge object --- ---

group <- factor(metadata$dicho_NIA_reagan) # Ajusta esto según tu metadata
dge <- DGEList(counts, group = group)

#Filter data on lower count rate 

drop <- filterByExpr(y = dge, min.count = 10, min.prop = 0.8)

dge <- dge[-drop,]
dim(dge)
#[1] 60557   793

#Normalize the data using the TMM scaling factor method.

dge <- calcNormFactors(dge)

# MDS plot
plotMDS(dge, labels = group, col = as.numeric(group))

# Boxplot de valores de cuentas normalizadas
logCPM <- cpm(dge, log = TRUE)
boxplot(logCPM, las = 2, col = as.numeric(group))

# Model Matrix 

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(dge, design)

# Differential expression 
colnames(design) <- c("dicho_NIA_Reagan_0", "dicho_NIA_Reagan_1")
rownames(design) <- metadata$specimenID
group1 <- "0"
group2 <- "1"
#voom for transformation and linear modeling
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)
#Use makeContrasts to define the contrasts you want to test
cont.matrix <- makeContrasts('dicho_NIA_Reagan_0 - dicho_NIA_Reagan_1', levels = design)
fit <- contrasts.fit(fit, cont.matrix )
fit <- eBayes(fit)
results <- topTable(fit, number = Inf, adjust.method = "BH")

#Vulcano plot --- ---

volc <- EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value')

#Heatmap --- ---

pheatmap(
  results,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Esquema de color (aquí azul-blanco-rojo)
  clustering_method = "complete",  # Método de agrupamiento jerárquico completo
  scale = "row",  # Escalar por filas (genes)
  show_rownames = TRUE,  # Mostrar nombres de filas (genes)
  show_colnames = TRUE,  # Mostrar nombres de columnas (condiciones/grupos)
  main = "Heatmap de Expresión Diferencial"  # Título del heatmap
)
