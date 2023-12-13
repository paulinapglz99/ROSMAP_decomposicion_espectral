#Script tidy-oso de preprocesamiento y QC de datos de expresion genica
#not finished

#paulinapglz.99@gmail.com

## &&
## Data preparation: 
##      -Quality Control & bias removal

#get the data

#Read normalized data (this input is for AD, MCI and noMCI, but here I use only AD)

dir <-'/datos/rosmap/FPKM_data/FPKM_AD.csv' #counts for every cogdx

FPKM <- vroom::vroom(file = dir)

#pull identifiers for annotation and delete version numbers

identifiers <- pull(FPKM, gene_id) %>% 
  substr(1, 15)

FPKM$gene_id <- identifiers

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

#annotate

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values= identifiers,
                 mart= mart)

#Add lenght column

myannot$length <- abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation

myannot <- myannot %>% 
  dplyr::filter(gene_biotype == "protein_coding" | hgnc_symbol!="") %>% #solo las filas donde ene_biotype sea "protein_coding" y hgnc_symbol no sea una cadena vacia 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Mantiene solo las filas únicas basadas en la columna ensembl_gene_id 

#Retener solo aquellas cuyos nombres de fila (rownames(expre)) coincidan con los valores deensembl_gene_id en el dataframe myannot.

FPKM_exprots <- FPKM %>%
  filter(FPKM$gene_id %in% myannot$ensembl_gene_id)
dim(FPKM_filtered)

#[1] 37963   256

#check duplicated probes

##################CHECK BIASES########################################################

pacman::p_load("NOISeq", 
               "edgeR")

#format data for noiseq

noiseqData <- readData(data = FPKM_exprots, 
                      gc = myannot[,1:2],           #porcentaje de GC en myannot
                      biotype = myannot[,c(1,3)],   #
                      factor = designExp,
                      length = myannot[,c(1,8)])


#1)check expression bias per subtype
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")
#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()
