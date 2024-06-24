#
#9.1 DEG_induced_subgraph.R
#This script filters a full edgelist to obtain nodes of the DEGS previously obtained

#Libraries --- ---
pacman::p_load("dplyr", 
               "vroom", 
               "igraph")

#Get data --- ---

#Differentially expressed genes
DEGS <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Graph
edgelist_0 <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_DEGS_dicho_edgelist.txt')

edgelist_1 <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_DEGS_dicho_edgelist.txt')

#Convert into graph 

graph_0 <- graph_from_data_frame(edgelist_0, directed = F)

graph_1 <- graph_from_data_frame(edgelist_1, directed = F)

#
graph_0.df <- data.frame(degree_g_0 = degree(graph_0))

graph_1.df <- data.frame(degree_g_1 = degree(graph_1))

degree <- cbind(graph_0.df,graph_1.df) %>% 
  rownames_to_column("ensembl_gene_id")

# Left join

DEGS <- DEGS %>% left_join(degree, by = "ensembl_gene_id")

#Plot scatter --- ---

#BaseMean

DEGS_BM_0 <- ggplot(DEGS, aes(x=degree_g_0, y=baseMean_NIA_R_0)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_0)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$baseMean_NIA_R_0)/2, col="red") +
  xlab("Degree") +
  ylab("baseMean") +
  labs(title = "", 
       subtitle = "NIA-Reagan 0")+
  theme_classic()
DEGS_BM_0

DEGS_BM_1 <- ggplot(DEGS, aes(x=degree_g_1, y=baseMean_NIA_R_1)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_1)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$baseMean_NIA_R_1)/2, col="red") +
  xlab("Degree") +
  ylab("baseMean") +
  labs(title = "", 
       subtitle = "NIA-Reagan 1") +
  theme_classic()
DEGS_BM_1

#logfoldChange

DEGS_lfc_0 <- ggplot(DEGS, aes(x=degree_g_0, y=log2FoldChange)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$log2FoldChange)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$log2FoldChange)/2, col="red") +
  xlab("Degree") +
  ylab("log2FoldChange") +
  labs(title = "", 
       subtitle = "NIA-Reagan 0") +
  theme_classic()

DEGS_lfc_1 <- ggplot(DEGS, aes(x=degree_g_1, y=log2FoldChange)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$log2FoldChange)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$log2FoldChange)/2, col="red") +
  xlab("Degree") +
  ylab("log2FoldChange") +
  labs(title = "", 
       subtitle = "NIA-Reagan 1") +
  theme_classic()

#Save plots



#END
