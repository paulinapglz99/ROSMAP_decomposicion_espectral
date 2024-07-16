#
# differential_expression_DEseq2_msex.R
#this script processes RNA-seq data 
#and performs the differential expression test to identify differentially expressed genes. 

#Differential expression between sick men and women
#Differential expression between healthy men and healthy women 

#Libraries --- ---

#BiocManager::install("DESeq2")

pacman::p_load("tidyverse", 
               "DESeq2", 
               "pheatmap", 
               "ggplot2", 
               "ggrepel", 
               "biomaRt", 
               "stringr",
               "RColorBrewer")

#Set seed --- --- 

set.seed(10)

#Define functions --- ---

#Function to translate gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connecting to the Ensembl database through biomaRt

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
}

#Get data --- ---

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
dim(counts)
#[1] 60603   880

#
rownames(counts) <- counts$feature
dim(counts)
#[1] 60603   1142

metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   41

table(metadata$dicho_NIA_reagan, useNA = "ifany")

#  0    1 <NA> 
#  307  573  261  

#Filter to obtain only the ones with NIA-Reagan dicho

metadata <- metadata %>% filter(!is.na(dicho_NIA_reagan))

#Only samples with metadata
counts <- counts %>% dplyr::select(all_of(metadata$specimenID))
dim(counts)
#[1] 60603   880

#Differential expression --- ---

#Experimental design

coldata <- as.data.frame(colnames(counts))
colnames(coldata) <- "specimenID"

coldata <- coldata %>% left_join(metadata, by = "specimenID")
coldata <- coldata %>% dplyr::select("specimenID", "sequencingBatch","msex","cogdx", "ceradsc", "dicho_NIA_reagan")
coldata <- coldata %>% mutate(condition = paste0("msex_",msex,"_NR_",dicho_NIA_reagan))
coldata$condition <- as.factor(coldata$condition)
coldata$dicho_NIA_reagan <- as.factor(coldata$dicho_NIA_reagan)  #Convert to factor

#DESeqData object --- ---

#3. Diff entre msex
##a. Filtrar conteos para tener unicamente enfermos en una y sanos en otra
##b. La comparacion seria entre hombres y mujeres enfermos & hombres y mujeres sanos

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ sequencingBatch + condition) #Comparar diagnostico y sexo

#Pre-filtering

# A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
smallestGroupSize <- 3
#Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
#The count of 10 is a reasonable choice for bulk RNA-seq.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)
#[1] 30657   880

#Specify conditions to compare
#Establishing 
dds$condition <- relevel(dds$condition, ref = "msex_0_NR_0")

dds <- DESeq(dds)  #Slow

#Results

#Contrast "msex_0_NR_0", "msex_0_NR1" (Healthy and sick women)

res_msex_0_NR_0_NR_1 <- results(dds, contrast = c("condition", "msex_0_NR_1", "msex_0_NR_0"), 
               pAdjustMethod = 'BH', 
               alpha = 0.01)

#Contrast "msex_1_NR_0", "msex_1_NR1" (Healthy and sick men)

res_msex_1_NR_0_NR_1 <- results(dds, contrast = c("condition", "msex_0_NR_1", "msex_0_NR_0"), 
                                pAdjustMethod = 'BH', 
                                alpha = 0.01)

#Contrast "msex_0_NR_1", "msex_1_NR_1"(Sick women and sick men)

res_NR_1_msex_0_msex_1 <- results(dds, contrast = c("condition", "msex_0_NR_1", "msex_1_NR_1"), 
               pAdjustMethod = 'BH', 
               alpha = 0.01)

#Contrast "msex_0_NR_1", "msex_1_NR_1"(Healthy women and Healthy men)

res_NR_1_msex_0_msex_1 <- results(dds, contrast = c("condition", "msex_0_NR_0", "msex_1_NR_0"), 
                                  pAdjustMethod = 'BH', 
                                  alpha = 0.01)

#Create a list 

res <- list(res_msex_0_NR_0_NR_1 = res_msex_0_NR_0_NR_1, 
            res_msex_1_NR_0_NR_1 = res_msex_1_NR_0_NR_1, 
            res_NR_1_msex_0_msex_1 = res_NR_1_msex_0_msex_1, 
            res_NR_1_msex_0_msex_1 = res_NR_1_msex_0_msex_1)

# Contar el número de p-valores ajustados menores a 0.05 en cada resultado
adjp <- sapply(res, function(x) sum(x$padj < 0.05, na.rm = TRUE))

res.f <- function(res) {
  #ordenar el data frame por padj
  res.df <- as.data.frame(res)
  res.df <- res.df[order(res.df$padj, decreasing = FALSE), ]
  # Convertir a data frame
  res.df <- as.data.frame(res)
  # Ordenar por log2FoldChange en orden decreciente
  res.df <- res.df[order(abs(res.df$log2FoldChange), decreasing = TRUE),]
  # Añadir una columna de NAs
  res.df$diffexpressed <- "NO"
  # Clasificar según log2FoldChange y p-valor ajustado
  res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
  res.df$diffexpressed[res.df$log2FoldChange < -0.5 & res.df$padj < 0.05] <- "DOWN"
  rownames(res.df) <- str_remove(rownames(res.df), "\\..*$")
  #Create dictionary
  symbol <- convert_ens_to_symbol(rownames(res.df))
  symbol$external_gene_name <- ifelse(symbol$external_gene_name == "", symbol$ensembl_gene_id, symbol$external_gene_name)
  #merge
  
  # Make rownames columns
  res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
  res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge
  #Create labels for Vplot
  res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA ))
  return(res.df)
}

# Función para extraer DEGs
DEGs.f <- function(res.df) {
  DEGS <- res.df %>% filter(diffexpressed != "NO")
  rownames(DEGS) <- DEGS$ensembl_gene_id
  return(DEGS)
}

# Aplicar la función a cada elemento de la lista
res.df <- lapply(res, res.f)

# Aplicar la función a cada elemento de la lista
DEGs.df <- lapply(res.df, DEGs.f)

# Mostrar las dimensiones de cada data frame de DEGs
DEGs_dims <- lapply(DEGs.df, dim)
DEGs_dims

# BaseMean by condition --- ---

#Extract counts normalized 

normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()

# Calcula los valores de baseMean por condición

# Definir las condiciones únicas
conditions <- unique(coldata$condition)

# Crear una lista para almacenar los data frames de baseMean
baseMeans_list <- list()

# Calcular los valores de baseMean para cada condición
for (condition in conditions) {
  baseMean <- rowMeans(normcounts[, coldata$condition == condition]) %>% as.data.frame()
  colname <- paste0("baseMean_", condition)
  colnames(baseMean) <- colname
  baseMeans_list[[colname]] <- baseMean
}

# Combinar los resultados en un solo data frame
baseMeans <- do.call(cbind, baseMeans_list)

# Add to every 

res.df <- lapply(res.df, function(df) {
  combined_df <- cbind(df, baseMeans[rownames(df), ])
  return(combined_df) })

#Vulcano plots

vplot.df <- 

#Create labels for Vplot
res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA ))

#
vplot <-  ggplot(data=a, aes(x=res_msex_0_NR_0_NR_1.log2FoldChange, y= -log10(res_msex_0_NR_0_NR_1.padj),  col = res_msex_0_NR_0_NR_1.diffexpressed, label = res_msex_0_NR_0_NR_1.delabel)) +
  geom_point() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(0.5, -0.5), col="red") + #log2FoldChange threshold is 0.5
  geom_hline(yintercept=-log10(0.05), col="red") + #p-value threshold is 0.05
  scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
  theme_minimal() +
  geom_text_repel() +
  labs(title = "Differential expression", 
       subtitle = "Using dichotomic NIA-Reagan Criteria")

b<- res.df[3] %>% as.data.frame()

vplot <-  ggplot(data=b, aes(x=res_NR_1_msex_0_msex_1.log2FoldChange, y= -log10(res_NR_1_msex_0_msex_1.padj),  col = res_NR_1_msex_0_msex_1.diffexpressed, label = res_NR_1_msex_0_msex_1.delabel)) +
  geom_point() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(0.5, -0.5), col="red") + #log2FoldChange threshold is 0.5
  geom_hline(yintercept=-log10(0.05), col="red") + #p-value threshold is 0.05
  scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
  theme_minimal() +
  geom_text_repel() +
  labs(title = "Differential expression", 
       subtitle = "Using dichotomic NIA-Reagan Criteria")

vplot

#Save log2fold change full list --- ---

#vroom::vroom_write(res.df, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_differential_expr_dichoNIAReagan.txt")

#Save DEGs list --- ---

#vroom::vroom_write(DEGS, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Create matrix of DEG counts --- ---

DEG_mat <- BiocGenerics::counts(dds, normalized = F) %>% as.data.frame()
rownames(DEG_mat)  <- str_remove(rownames(DEG_mat) , "\\..*$")

DEG_mat <- DEG_mat %>% filter(rownames(DEG_mat) %in% DEGS$ensembl_gene_id)
dim(DEG_mat)
#[1] 4130  880

# Create correlation matrix (distance matrix)

# DEG_mat.z <- t(apply(DEG_mat, 1, scale))
# #Add sample names
# colnames(DEG_mat.z) <- colnames(DEG_mat) 
# 
# DEG_mat.z[1:20, 1:20]

#You can also try (I think this is better)

#Estimate dispersion trend and apply a variance stabilizing transformation
vsd <- vst(DEG_mat, 
           blind = F) #whether to blind the transformation to the experimental design 
vsd_assay <- assay(vsd) %>% t() #assay function is used to extract the matrix of normalized values
colnames(vsd_assay) <- str_remove(colnames(vsd_assay) , "\\..*$")

#Sample distance metrix
sampleDists <- dist(vsd_assay) #slow
sampleDistsMatrix <- as.matrix(sampleDists) #

rownames(sampleDistsMatrix) <- paste(vsd$specimenID, vsd$dicho_NIA_reagan, sep="-")
colnames(sampleDistsMatrix) <- paste(vsd$specimenID, vsd$dicho_NIA_reagan, sep="-")

# Vis

sampleDistsMatrix[1:5, 1:5]

#Heatmap --- ---

# Extract sample information
sample_info <- data.frame(
  AD_NIA_Reagan_dichotomous = vsd$dicho_NIA_reagan
)
rownames(sample_info) <- rownames(sampleDistsMatrix)

condition_colors <- list(
  AD_NIA_Reagan_dichotomous = c("0" = "blue", "1" = "red") # adjust the colors and names according to your data
)

#Heatmap of sample-to-sample distance matrix (with clustering) based on normalized counts ---

#Set a color scheme
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)

#Generate the heatmap

#This heatmap tells you what samples are more likely to each other

distances.p <- pheatmap(sampleDistsMatrix,
                        annotation_col = sample_info, 
                        annotation_row = sample_info, 
                        annotation_colors = condition_colors, 
                        show_rownames = F, 
                        show_colnames = F, 
                        clustering_distance_rows=sampleDists,
                        clustering_distance_cols=sampleDists,
                        col=colors)

distances.p

#Heatmap of log transformed normalized counts with top10 genes ---

# Choose top 20 genes
DEG_top <- res.df %>% filter(!is.na(padj))
DEG_top <- DEG_top[order(res.df$padj),][1:20,]
rownames(DEG_top) <- DEG_top$ensembl_gene_id

#Log transformation

rld <- rlog(dds, blind = FALSE) #slow

annot_info <- as.data.frame(colData(dds)[, c("msex", "dicho_NIA_reagan")])
dim(annot_info)

tr<- vsd_assay[ ,DEG_top10$ensembl_gene_id]
dim(tr)

#plot

DEG_top10.p <- pheatmap(tr,
                        cluster_rows = F, 
                        show_rownames = T, 
                        cluster_cols = F, 
                        annotation_col = annot_info, 
                        annotation_row = annot_info)

#Heatmap of Z scores. Top 10 genes --- 

# Get normalized counts

norm_counts <- counts(dds, normalized=TRUE)

#The Z-score gives the number of standard-deviations that a value is away from the mean of all the values in the same group, here the same gene. 

#Compute the Z-scores 

zscore_all <- t(apply(norm_counts, 1, cal_z_score))

zscore_subset <- zscore_all[DEG_top, ]

#Plot

zscores.p <- pheatmap(zscore_all)

#with complexheatmap ---

col_logFC <- colorRamp2(c(min()))

ComplexHeatmap::Heatmap(DEG_mat.z, cluster_rows = T, cluster_columns = T, 
                        column_labels = colnames(DEG_mat.z),
                        column_title = "Samples", row_title = "DEGs",
                        # column_split = list(specimenID = condition$specimenID, dicho_NIA_reagan = condition$dicho_NIA_reagan),
                        name = "Z-score", row_labels = DEGS[rownames(DEG_mat.z),]$external_gene_name)

# Save heatmaps --- ---

#distances.p
ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_sample_distances.png", 
       plot = distances.p, 
       width = 11,
       height = 15,
       units = "in",
       dpi = 300)

#DEG_top10.p

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_DEG_top10.png", 
       plot = DEG_top10.p, 
       width = 11,
       height = 15,
       units = "in",
       dpi = 300)

#zscores.p

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_zscores.png", 
       plot = zscores.p, 
       width = 11,
       height = 15,
       units = "in",
       dpi = 300)

#Vulcano plot --- ---

#Create labels for Vplot
res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA))

#
vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y=-log10(pvalue),  col = diffexpressed, label = delabel)) +
  geom_point() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") + #log2FoldChange threshold is 0.6
  geom_hline(yintercept=-log10(0.05), col="red") + #p-value threshold is 0.05
  scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
  theme_minimal() +
  geom_text_repel() +
  labs(title = "Differential expression", 
       sub = "Using dichotomic NIA-Reagan Criteria")

#Save Vulcano Plot ---

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_vulcano_plot.png", 
       plot = vplot, 
       width = 11,
       height = 15,
       units = "in",
       dpi = 300)
