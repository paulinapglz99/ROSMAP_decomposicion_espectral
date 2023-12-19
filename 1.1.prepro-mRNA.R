#Script tidy preprocessing and QC of gene expression data

#paulinapglz.99@gmail.com

## Data preparation: 
## Quality Control & bias removal

pacman::p_load('dplyr', 
               'biomaRt')

############################## A. Get the data ############################## 

#Read normalized data (this input is for AD, MCI and noMCI, but here I use only AD)
#This file was generated in 1.MatchFPKMandClinicalMetadata.R

FPKM <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') #counts for cogdx = 2, 3, 4 and 5
dim(FPKM)
#[1] 55889   623 #for AD, original expression counts have 55889 genes and 623 specimenIDs

############################## B. Annotation ##############################

#pull identifiers for annotation and delete version numbers

identifiers <- pull(FPKM, gene_id)

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

#annnotate GC content, length & biotype per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values = identifiers,
                 mart = mart)

#Add lenght column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49410     8   #The full annotation has 49410 genes

#filter transcripts withouth annotation

myannot <- myannot %>% 
  dplyr::filter(gene_biotype == "protein_coding" | hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column
dim(myannot)
#[1] 37963     8  #37963 full annotated protein coding genes with 8 characteristics, difference of 11447 annotated genes :c

#Filter count matrix to have only values that match ensembl_gene_id in the myannot dataframe.

FPKM_exprots <- FPKM %>%
  filter(FPKM$gene_id %in% myannot$ensembl_gene_id)
dim(FPKM_exprots)
#[1] 37963   623  #37963 genes from 622 specimenIDs

################## C. CHECK BIASES ########################################################

pacman::p_load("NOISeq", 
               "edgeR")

#assign design experiment format

#read metadata
RNA_seq_metadata <- vroom::vroom(file = "/datos/rosmap/metadata/ROSMAP_filtered_metadata_forRNAseq.csv")
dim(RNA_seq_metadata)
#[1] 622   5   #622 individualIDs and specimenIDs and 5 characteristics

#The lengths of RNA_seq_metadata$specimenID and colnames(FPKM_exprots) must be the same for the noiseq::readData function

# Extraer todas las columnas excepto la primera de FPKM_exprots FOR NOISEQ
common_samples <- intersect(colnames(FPKM_exprots), RNA_seq_metadata$specimenID)
#chr [1:622]

#here we see that both lenghts match, (omitting the gene_id column of course)

FPKM_exprots <-FPKM_exprots[, common_samples]
dim(FPKM_exprots)
#[1] 37963   622 ## is the same in this case

RNA_seq_metadata <- RNA_seq_metadata %>% 
  filter(specimenID %in% colnames(FPKM_exprots))
dim(RNA_seq_metadata)
#[1] 622   5 #esto ya coincide

#for the factor format
#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data le provided.

# Get indexes to reorder RNA_seq_metadata$specimenID
specimenID_positions <- match(colnames(FPKM_exprots), RNA_seq_metadata$specimenID)

# Use indexes to reorder RNA_seq_metadata$specimenID
RNA_seq_metadata$specimenID <- RNA_seq_metadata$specimenID[specimenID_positions]
names(RNA_seq_metadata$specimenID) <- colnames(FPKM_exprots)

# Now, RNA_seq_metadata$specimenID_ordered contains the values of RNA_seq_metadata$specimenID
# sorted according to the order of colnames(FPKM_exprots)

#note: row number in factor must match with number of cols in data ("FPKM_exprots"). 

noiseqData <- readData(data = FPKM_exprots, 
                      gc = myannot[,c("ensembl_gene_id","percentage_gene_gc_content")],  #%GC in myannot
                      biotype = myannot[,c("ensembl_gene_id", "gene_biotype")],          #biotype
                      factor = RNA_seq_metadata,                                         #variables indicating the experimental group for each sample
                      length = myannot[,c("ensembl_gene_id", "length")])                 #gene length

#1)check expression bias per subtype

mycountsbio = dat(noiseqData,
                  type = "countsbio",
                  factor = "subtype")

#patients with repeated measures

png("CountsOri.png")

explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)

##Next script is 1.2.
