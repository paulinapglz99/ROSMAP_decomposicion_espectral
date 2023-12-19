#Script tidy-oso de preprocesamiento y QC de datos de expresion genica

#paulinapglz.99@gmail.com

## Data preparation: 
## Quality Control & bias removal

############################## A. Get the data ############################## 

#Read normalized data (this input is for AD, MCI and noMCI, but here I use only AD)

FPKM <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') #counts for cogdx = 2, 3, 4 and 5

#pull identifiers for annotation and delete version numbers

identifiers <- pull(FPKM, gene_id)

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

#annotate

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values = identifiers,
                 mart = mart)

#Add lenght column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49410     8

#filter transcripts withouth annotation

myannot <- myannot %>% 
  dplyr::filter(gene_biotype == "protein_coding" | hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column
dim(myannot)
#[1] 37963     8  #37963 genes in 8 annotations

#Retener solo aquellas cuyos nombres de fila (rownames(expre)) coincidan con los valores deensembl_gene_id en el dataframe myannot.

FPKM_exprots <- FPKM %>%
  filter(FPKM$gene_id %in% myannot$ensembl_gene_id)
dim(FPKM_exprots)

#[1] 37963   256

#check duplicated probes

##################CHECK BIASES PAU ########################################################

pacman::p_load("NOISeq", 
               "edgeR")

#assign design experiment

RNA_seq_metadata <- vroom::vroom(file = "/datos/rosmap/metadata/ROSMAP_filtered_metadata_designExp.csv")
dim(designExp)

format <- RNA_seq_metadata %>%  #necesito hacer que la longitud de la columnna specimenID de 
          select()                   #esta tabla coincida con el numero (y nombre) de columas de mi matriz de conteos

#filas columnas
#[1] 634   5

#format data for noiseq
#note: row number in factor ("designExp") must match with number of cols in data ("FPKM_exprots"). 

noiseqData <- readData(data = FPKM_exprots, 
                      gc = myannot[,c("ensembl_gene_id","percentage_gene_gc_content")],           #porcentaje de GC en myannot
                      biotype = myannot[,c("ensembl_gene_id", "gene_biotype")],   #
                      factor = designExp, 
                      length = myannot[,c("ensembl_gene_id", "length")])

#format data for noiseq

noiseqData = readData(data = exprots_hgnc, gc = myannot[,1:2],
                      biotype = myannot[,c(1,3)],
                      factor = designExp,
                      length=myannot[,c(1,8)])


#1)check expression bias per subtype
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")
#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()


##Next script is 1.2.
