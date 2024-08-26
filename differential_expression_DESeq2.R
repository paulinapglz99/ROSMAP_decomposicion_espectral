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

counts <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_filteredQC_counts_DLPFC.rds") %>% as.data.frame()
dim(counts)

# Get metadata
metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt")
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
coldata <- coldata %>% dplyr::select("specimenID", "sequencingBatch","sex","cogdx", "ceradsc", "dicho_NIA_reagan")
coldata$dicho_NIA_reagan <- as.factor(coldata$dicho_NIA_reagan)  #Convert to factor

#DESeqData object

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = coldata,
                              design = ~ dicho_NIA_reagan) #En caso unicamente de comparar diagnostico

#Differential expression analysis --- ---

#Specify conditions to compare
#Establishing 0 (not AD) as reference
dds$dicho_NIA_reagan <- relevel(dds$dicho_NIA_reagan, ref = "0")

dds <- DESeq(dds)  #Slow

#Results

res <- results(dds,
               contrast = c("dicho_NIA_reagan", "1", "0"), #contrast = c("condition", "problem", "control"))
               pAdjustMethod = 'BH',   #he method to use for adjusting p-values
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
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
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
#vroom::vroom_write(DEGS, file = "q")

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
  geom_text_repel(max.overlaps = 50) +
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

DEG_mat <- BiocGenerics::counts(dds, normalized = T) %>% as.data.frame()
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

normcounts_DEG <- normcounts[rownames(DEG_mat),]

htmap <- pheatmap(DEG_mat,
         cluster_rows = T,   # Cluster the rows
         cluster_cols = T,   # Do not cluster the columns
         main = "Heatmap of Differentially Expressed Genes", 
     #    annotation_col = metadata_DEG,
         show_colnames = FALSE  # Hide the column names (sample names)
)

# Save heatmaps --- ---
 
# ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_zscores.png", 
#        plot = htmap, 
#        width = 11,
#        height = 15,
#        units = "in",Screenshot from 2024-08-21 15-17-00
#        dpi = 300)