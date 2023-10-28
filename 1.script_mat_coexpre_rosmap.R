#This script takes FPKM normalized data from RNAseq and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block

#Libraries

pacman::p_load("tidyverse", 
       "ggplot2", 
       'vroom', 
       'biomaRt')

#Read data

dir<-'FPKM_AD.csv' #discretized coexpression matrix

FPKM <- vroom(file = dir)

##join tables with different expression data plates IF NEEDED

#FPKM <-dplyr::left_join(x=genes_expre_p1_p6,
#                              y=genes_expre_p7_p8,
#                             by=c('gene_id'='gene_id'))
#
#FPKM$tracking_id.x=NULL

#vroom_write(FPKM, 
#           file = 'RNAseq_FPKM_1_to_8_merged.csv', 
#           delim = ',')

##Clean column 'gene_id': remove .XX from gene_id

identificadores<- pull(FPKM,gene_id)

identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])

FPKM <- FPKM %>% add_column(identificadores)

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", 
                                "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values= identificadores,
                 mart= mart)

##join biomart annotations with data table, to generate the annotated data

expre <-left_join(x = FPKM,
                  y = myannot,
                  by = c('identificadores'='ensembl_gene_id'))

##Keep only the 'gene coding' type data (since for now they are the only ones we are interested in).

protcod <- filter(expre, 
                       gene_biotype =='protein_coding')

valores_expre <- protcod %>% 
  dplyr::select(-gene_id,
                -gene_biotype,
                -identificadores,
                -percentage_gene_gc_content)

#### Discretization of data
mat_dis<-infotheo::discretize(t(valores_expre))

#Save matrix

#vroom_write(mat_dis, 
#       file = 'coexpression_matrix_AD.txt', 
#      delim = ',')

###Next script is 3.mutualinformation_parallel.R