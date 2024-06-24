#
#3.differential_expression_DEseq2.R
#this script processes RNA-seq data to identify differentially expressed genes. 

#Libraries --- ---

#BiocManager::install("DESeq2")

pacman::p_load("tidyverse", 
               "DESeq2", 
               "pheatmap", 
               "ggplot2", 
               "ggrepel", 
               "stringr",
               "biomaRt", 
               "RColorBrewer")

#Set seed --- --- 

set.seed(10)

#Define functions --- ---

#Function to translate gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connecting to the Ensembl database through biomaRt

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
 trad <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
  trad$external_gene_name <- ifelse(trad$external_gene_name == "", trad$ensembl_gene_id, trad$external_gene_name)
  return(trad)
}

#Get data --- ---

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
counts$feature <- str_remove(counts$feature, "\\..*$")

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
#  Convert selected columns to integers using dplyr
counts <- counts %>% mutate(across(-feature, as.integer))

# Get metadata
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
#[1] 60558   880

#Differential expression --- ---

#Experimental design

coldata <- as.data.frame(colnames(counts))
colnames(coldata) <- "specimenID"

coldata <- coldata %>% left_join(metadata, by = "specimenID")
coldata <- coldata %>% dplyr::select("specimenID", "sequencingBatch","msex","cogdx", "ceradsc", "dicho_NIA_reagan")
coldata$dicho_NIA_reagan <- as.factor(coldata$dicho_NIA_reagan)  #Convert to factor

#DESeqData object

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = coldata,
                              design = ~ sequencingBatch + dicho_NIA_reagan) #En caso unicamente de comparar diagnostico

#Pre-filtering

# A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
smallestGroupSize <- 3
#Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
#The count of 10 is a reasonable choice for bulk RNA-seq.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)
#[1] 30657   880

#Differential expression analysis --- ---

#Specify conditions to compare
#Establishing 0 (not AD) as reference
dds$dicho_NIA_reagan <- relevel(dds$dicho_NIA_reagan, ref = "0")

dds <- DESeq(dds)  #Slow

#Results

res <- results(dds,
               contrast = c("dicho_NIA_reagan", "1", "0"), #contrast = c("condition", "problem", "control"))
               pAdjustMethod = 'BH', 
               alpha = 0.05) #FDR

#How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)
#[1] 4423

res.df <- res %>% as.data.frame()
res.df <- res.df %>% filter(!is.na(padj))
dim(res.df)
#[1] 27683     6

# add a column of NAs
res.df$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.05, set as "DOWN"
res.df$diffexpressed[res.df$log2FoldChange < -0.5 & res.df$padj < 0.05] <- "DOWN"

table(res.df$diffexpressed)

#DLPFC log2FoldChange > 0.5 & res.df$padj < 0.05
# 
# DOWN    NO    UP 
#   8 27655    20 

#Add gene names in SYMBOL --- ---

#Create dictionary
symbol <- convert_ens_to_symbol(rownames(res.df))

# Make rownames columns
res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge

# BaseMean by condition --- ---

#Extract counts normalized 

normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()

# Calcula los valores de baseMean por condición
baseMean_0 <- rowMeans(normcounts[, coldata$dicho_NIA_reagan == "0"]) %>% as.data.frame()
colnames(baseMean_0) <- "baseMean_NIA_R_0"

baseMean_1 <- rowMeans(normcounts[, coldata$dicho_NIA_reagan == "1"])%>% as.data.frame()
colnames(baseMean_1) <- "baseMean_NIA_R_1"

baseMean_genes <- cbind(baseMean_0, baseMean_1)
baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")

#Extract DEGs --- ---

DEGS <- res.df %>% filter(diffexpressed != "NO")
rownames(DEGS) <- DEGS$ensembl_gene_id
dim(DEGS)
#[1] 28    10

baseMean_DEGS <- baseMean_genes %>% filter(feature %in%rownames(DEGS))
DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))

#Save log2fold change full list --- ---
# 
#vroom::vroom_write(res.df, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_differential_expr_dichoNIAReagan.txt")
# 
# #Save DEGs list --- ---
# 
#vroom::vroom_write(DEGS, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Vulcano plot --- ---

#Create labels for Vplot
res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA ))

#
vplot <-  ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj),  col = diffexpressed, label = delabel)) +
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
#Save Vulcano Plot ---

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_vulcano_plot.png", 
       plot = vplot, 
       width = 11,
       height = 15,
       units = "in",
       dpi = 300)

#Create matrix of DEG counts --- ---

DEG_mat <- BiocGenerics::counts(dds, normalized = F) %>% as.data.frame()
DEG_mat <- DEG_mat %>% filter(rownames(DEG_mat) %in% DEGS$ensembl_gene_id)
dim(DEG_mat)
#[1] 28  880

# Create correlation matrix (distance matrix)

# DEG_mat.z <- t(apply(DEG_mat, 1, scale))
# #Add sample names
# colnames(DEG_mat.z) <- colnames(DEG_mat) 
# 
# DEG_mat.z[1:20, 1:20]

#You can also try (I think this is better)

#Estimate dispersion trend and apply a variance stabilizing transformation
rld <- rlog(as.matrix(DEG_mat), 
           blind = F) #whether to blind the transformation to the experimental design 

trad <- convert_ens_to_symbol(rownames(rld))

rownames(rld) <- trad$external_gene_name

# Z score heatmap --- ---

Z <- t(scale(t(rld)))
dim(Z)

# Crear el heatmap

z.p <- pheatmap(Z, 
         name = "Row Z-Score", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         clustering_distance_rows = "euclidean", 
         clustering_method = "ward.D2", 
         fontsize = 11, 
         angle_col = "45", 
         show_rownames = T)

# Save heatmaps --- ---

#zscores.p
# 
# ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_zscores.png", 
#        plot = zscores.p, 
#        width = 11,
#        height = 15,
#        units = "in",
#        dpi = 300)