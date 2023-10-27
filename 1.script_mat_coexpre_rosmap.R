#Este script toma datos normalizados en FPKM de RNAseq y calcula una matriz de coexpresion
#Con esta matriz de coexpresion se realiza mas adelante la matriz de informacion mutua

#paquetes

pacman::p_load("tidyverse", 
       "ggplot2", 
       'vroom', 
       'biomaRt')

##abrir archivos
dir<-'ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv'  #archivo
dir2 <-'ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv'
genes_expre_p1_p6 <- vroom(file = dir)
genes_expre_p7_p8 <- vroom(file = dir2)

##join tablas con datos de expresión diferentes plates SI ES NECESARIO

expre_junto<-dplyr::left_join(x=genes_expre_p1_p6,
                              y=genes_expre_p7_p8,
                              by=c('gene_id'='gene_id'))

expre_junto$tracking_id.x=NULL


#vroom_write(expre_junto, 
#           file = 'RNAseq_FPKM_1_to_8_merged.csv', 
#           delim = ',')

##Limpiar columna 'gene_id': quitar los puntos del gene_id (tomo parte del código de sol)

identificadores<-expre_junto %>% 
  pull(gene_id) 

identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])

expre_junto<-expre_junto %>% add_column(identificadores)

#Generar anotacion con ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

myannot <- getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", 
                             "gene_biotype"),
              filters = "ensembl_gene_id", 
              values= identificadores,
              mart= mart)

##juntar anotaciones de biomart con tabla de datos, para generar la data anotada

expre <-left_join(x = expre_junto,
                         y = myannot,
                         by = c('identificadores'='ensembl_gene_id'))

##Quedarme solo con los datos de tipo 'gene coding' (pues por ahora son los únicos que nos interesan)

solo_protcod <- filter(expre, 
                       gene_biotype =='protein_coding')

valores_expre <-solo_protcod %>% 
  dplyr::select(-gene_id,
                -gene_biotype,
                -identificadores,
                -percentage_gene_gc_content)



####  Discretizacion de los datos
mat_dis<-infotheo::discretize(t(valores_expre))

###De aqui el siguiente sript es script_calculo_mi_parallel.R