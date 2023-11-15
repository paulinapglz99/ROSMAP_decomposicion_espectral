#This script takes FPKM normalized data from RNAseq and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block

#Libraries

pacman::p_load("tidyverse", 
       "ggplot2", 
       'vroom', 
       'biomaRt')

#Read normalized data (this input is for AD, MCI and noMCI, but here I use only AD

dir <-'FPKM_AD.csv' #counts for every cogdx

FPKM <- vroom(file = dir)

#pull identifiers for annotation

identifiers <- pull(FPKM,gene_id)

#

identifiers <-sapply(strsplit(identifiers,".",fixed=T),function(x) x[1])

FPKM <- FPKM %>% add_column(identifiers)

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", 
                                "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values= identifiers,
                 mart= mart)

##join biomart annotations with data table, to generate the annotated data

expre <-left_join(x = FPKM,
                  y = myannot,
                  by = c('identifiers'='ensembl_gene_id'))

##Keep only the 'gene coding' type data (since for now they are the only ones we are interested in).

protcod <- filter(expre, 
                       gene_biotype =='protein_coding')

valores_expre <- protcod %>% 
  dplyr::select(-gene_id,
                -gene_biotype,
                -identifiers,
                -percentage_gene_gc_content)

#### Discretization of data
mat_dis<-infotheo::discretize(t(valores_expre))

#Save matrix

#vroom_write(mat_dis, 
#       file = 'coexpression_matrix_AD.txt', 
#      delim = ',')

###Next script is 3.mutualinformation_parallel.R
