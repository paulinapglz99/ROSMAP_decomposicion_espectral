#This script takes FPKM normalized data from RNAseq and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block
#Previous script is 1.1.pre-promRNA.R

#Libraries

pacman::p_load("tidyverse", 
       "ggplot2", 
       'vroom', 
       'biomaRt')


#Subset by cognitive diagnosis ------------------

#we have two diagnosis indicators in our metadata: cogdx and ceradsc
#cogdx is blind to post morten data
#ceradsc is blind to clinical data (histopatology)

#as the metadata dictionary says:

#for cogdx 
##1 NCI: No cognitive impairment (No impaired domains)
##2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
##3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
##4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD) 
##5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
##6 Other dementia: Other primary cause of dementia

#for ceradsc

#value | coding | if using a binary variable, recommendation is
#  1  | Definite | yes
#  2  | Probable | yes
#  3  | Possible | no
#  4  |  No AD   | no

#only if we have ceradsc = 1 or 2 AND cogdx = 4 we can say we have a fully confirmed AD diagnosis

#First, we assign libraries for the respective cognitive diagnosis

#In 1.MatchFPKMandClinicalMetadata i already deleted cogdx == 1 and 6

RNA_seq_metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_filtered_metadata_forRNAseq.csv')

#library for no MCI (no dementia known at the moment of death)

noMCIcogdx <- RNA_seq_metadata %>% 
  filter(cogdx == 1)
dim(noMCIcogdx)
# [1] 200   5   #200 individualIDs with no dementia

# For all MCI (mild cognitive impairment)

allMCIcogdx <- RNA_seq_metadata %>% 
  filter(cogdx == 2 | cogdx == 3)
dim(allMCIcogdx)
#[1] 168   5    #168 individuals with some type of MCI

# Library for all types of AD

allADcogdx <- RNA_seq_metadata %>% 
  filter(cogdx == 4 | cogdx == 5)
dim(allADcogdx)
#[1] 254   5    #254 individuals with probable or possible AD

#Alternatively, if we want only fully confirmed AD patients
#we need to choose by ceradsc == 1 OR 2 && cogdx == 4

confirmedAD <- RNA_seq_metadata %>% 
  filter(ceradsc %in% c(1, 2), cogdx == 4)
dim(confirmedAD)
#[1] 192   5    #192 invididuals with confirmed AD

#Difference of 62 individualIDs between allcogdxAD and confirmedAD

#Read count data, this is now QC 

#############################WARNING THIS IS NOT THE FINAL FILE, I HAVENT DONE THE QC YET ##########
fpkm_matrix <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') ###
#############################WARNING THIS IS NOT THE FINAL FILE, I HAVENT DONE THE QC YET ##########

# Finally we subset RNAseq FPKMS by cogdx ----------

FPKM_noMCI <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(noMCIcogdx$specimenID))
dim(FPKM_noMCI)
#[1] 55889   201  #201 specimenIDs

FPKM_MCI <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(allMCIcogdx$specimenID))
dim(FPKM_MCI)
# [1] 55889   169  #169 specimenIDs

FPKM_AD <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(allADcogdx$specimenID))
dim(FPKM_AD)
#[1] 55889   255  #255 spcimenIDs

FPKM_confirmedAD <- fpkm_matrix %>% 
  dplyr::select(gene_id, all_of(confirmedAD$specimenID))
dim(FPKM_confirmedAD)
#[1] 55889   193  #193 specimenIDs

#Save tables

#vroom_write(FPKM_AD, 
#        file = 'FPKM_AD.csv', 
#       delim = ',')

FPKM <- FPKM_AD  #For non-hardcoding reasons, change this to FPKM_MCI, FPKM_noMCI, FPKM_confirmedAD, etc.

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", 
                                "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values = FPKM$gene_id,
                 mart = mart)

##join biomart annotations with data table, to generate the annotated data

expre <-left_join(x = FPKM,
                  y = myannot,
                  by = c('gene_id'='ensembl_gene_id'))                

##Keep only the 'gene coding' type data (since for now they are the only ones we are interested in).

protcod <- filter(expre, 
                  gene_biotype =='protein_coding')
dim(protcod)
#[1] 19064   257

#delete name column, we dont need it for discretization

valores_expre <- protcod %>% 
  dplyr::select(-gene_id, 
                -percentage_gene_gc_content, 
                -gene_biotype)
dim(valores_expre)
#[1] 19064   254  #19064 genes

#pu

gene_names <- pull(protcod, 'gene_id')

vroom_write(protcod, 
            file = '/datos/rosmap/discretized_matrix/protcod_AD.txt',  #this for MCI and noMCI also
            delim = ',')

####  Data discretization
#this generates a discretized expression matrix

mat_dis<-infotheo::discretize(t(valores_expre))
dim(mat_dis)
#[1]   254 19064  #still 19064 genes yay

###Next script is 3.mutualinformation_matrix.R
