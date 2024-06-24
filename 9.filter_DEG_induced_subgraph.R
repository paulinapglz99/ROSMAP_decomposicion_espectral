#
#9.filter_DEG_induced_subgraph.R
#This script filters a full edgelist to obtain nodes of the DEGS previously obtained

#Libraries --- ---
pacman::p_load("dplyr", 
               "vroom")

#Get data --- ---

#Differentially expressed genes
DEGS <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Graph
edgelist <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.txt")        
colnames(edgelist)

#Filter edgelist

edgelist_filtered <- edgelist %>%
  filter(feature %in% DEGS$ensembl_gene_id | gene_to %in% DEGS$ensembl_gene_id )
dim(edgelist_filtered)

#Save DEG edgelist --- ---

vroom::vroom_write(edgelist_filtered, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_DEGS_dicho_edgelist.txt")

#END