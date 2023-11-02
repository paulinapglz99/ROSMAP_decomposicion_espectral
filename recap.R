#what ive done recap

#disclaimer: este codigo esta cochino, solo fue para poner en orden mi workflow

#pasos a partir de los datos crudos de rosmap

Sys.time()

#librerias 

pacman::p_load("tidyverse", 
               "ggplot2", 
               'vroom', 
               'biomaRt')

#lectura de datos

#metadata

clinical_metadata <- vroom(file = 'ROSMAP_clinical.csv')

biospecimen_metadata <- vroom(file = 'ROSMAP_biospecimen_metadata.csv')

cli_bio_metadata <- left_join(x = clinical_metadata, 
                              y = biospecimen_metadata, 
                              by = "individualID")


RNA_seq_metadata <- cli_bio_metadata %>% 
  filter(assay == 'rnaSeq',
         organ == 'brain')  #de una vez

#Match Clinical.csv file to RNA Seq files and subset 


#data

FPKM_p1_p6 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv')
FPKM_p7_p8 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv')

FPKM_p1_p8<-left_join(x=FPKM_p1_p6,
                      y=FPKM_p7_p8,
                      by=c('gene_id'='gene_id'))

# gene ID, the copy of the gene ID and two repeated samples with _6 and _7 are deleted (_0 was left)

FPKM_p1_p8 <- FPKM_p1_p8 %>% 
  dplyr::select(-c("tracking_id.x", "tracking_id.y", 
            "492_120515_6", "492_120515_7"))

#note: I lost 2 samples (492_120515_6 and 492_120515_7), and removed column tracking_id

#Remove the last two characters from geneID identifiers to match them with 
#the RNASeq metadata

colnames_woID <- substr(colnames(FPKM_p1_p8),1, nchar(colnames(FPKM_p1_p8))-2)
colnames_woID[1] <- "gene_id"  #as we remove the last two characters we have to recompose the gene_id

#subset por cogdx

fpkm_matrix <- FPKM_p1_p8[,(colnames_woID %in% RNA_seq_metadata$specimenID)] %>% 
  mutate(gene_id = FPKM_p1_p8$gene_id, .before = 1) %>%
  rename_at(-1, ~str_sub(., end = -3))   # Elimina los dos últimos caracteres de los nombres de columna menos de la gene_id

#

RNA_seq_metadata <- RNA_seq_metadata[(RNA_seq_metadata$specimenID %in% colnames_woID),c("specimenID","individualID","cogdx")] %>%  # 19 = specimenID 1 = individualID 17 = cogdx
  arrange(match(specimenID, colnames_woID[-1]))


#1 NCI: No cognitive impairment (No impaired domains)
#2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
#3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
#4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)
#5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
#6 Other dementia: Other primary cause of dementia

#we dont need cogdx == 6

#library for cogdx == 1 (NCI)

cogdx1 <- RNA_seq_metadata %>% 
  filter(cogdx == 1)

cogdx1v <- pull(cogdx1, specimenID)

# cogdx == c(2, 3) (all MCI)

cogdx2_3 <- RNA_seq_metadata %>% 
  filter(cogdx == 2 | cogdx == 3)

cogdx2_3v <- pull(cogdx2_3, specimenID)

# cogdx == c(4, 5) (all AD)

cogdx4_5 <- RNA_seq_metadata %>% 
  filter(cogdx == 4 | cogdx == 5)

cogdx4_5v <- pull(cogdx4_5, specimenID)


###########Finally we subset RNAseq FPKMS by cogdx

FPKM_noMCI <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(cogdx1v))

FPKM_MCI <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(cogdx2_3v))

FPKM_AD <- fpkm_matrix %>%
  dplyr::select(gene_id, all_of(cogdx4_5v))


#Save tables

#vroom_write(FPKM_AD, 
#        file = 'FPKM_AD.csv', 
#       delim = ',')


######### A partir de aqui se hace para cada tabla de FPKM


#filtrado para tener solo los protein coding

identificadores<- pull(FPKM_AD,gene_id)

identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])

FPKM_AD <- FPKM_AD %>% add_column(identificadores)

#Generar anotacion con ensembl. Annotate gene_biotype

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")


myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values = identificadores,
                 mart = mart)

##juntar anotaciones de biomart con tabla de datos, para generar la data anotada

expre <-left_join(x = FPKM_AD,
                  y = myannot,
                  by = c('identificadores'='ensembl_gene_id'))

##Quedarme solo con los datos de tipo 'gene coding' (pues por ahora son los únicos que nos interesan)

protcod <- filter(expre, 
                  gene_biotype =='protein_coding')

valores_expre <- protcod %>% 
  dplyr::select(-gene_id,
                -identificadores)

my_index <- pull(protcod, 'gene_id')

#vroom_write(protcod, 
#            file = 'protcod_AD.txt', 
#           delim = ',')

####  Discretizacion de los datos
#esto genera una matriz de expresion discretizada

mat_dis<-infotheo::discretize(t(valores_expre))


########################################

#A partir de aqui, es el script de Guillermo 

plan(multicore, workers = 40)

# calculate MI

mi_mi <- 
  map(my_index, .f = function(i){
  ii = mat_dis[i]
  map(my_index, .f = function(j){
    jj = mat_dis[j]
    #if(j>i)
    {
      
      mutinformation(ii,jj)
      
    }
    
  })
  
})

#write out ---- 

write_rds(x =mi_mi, file = "/datos/rosmap/MImatrix_AD.rds")

print(Sys.time() - tempus())
