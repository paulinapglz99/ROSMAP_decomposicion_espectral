#Este script toma datos normalizados en FPKM de RNAseq y calcula una matriz de coexpresion
#Con esta matriz de coexpresion se realiza mas adelante la matriz de informacion mutua

#paquetes

pacman::p_load("tidyverse", 
       "ggplot2", 
       'vroom', 
       'biomaRt')

##abrir archivos

dir<-'FPKM_AD.csv' #archivo

FPKM <- vroom(file = dir)

##join tablas con datos de expresión diferentes plates SI ES NECESARIO

#expre_junto<-dplyr::left_join(x=genes_expre_p1_p6,
#                              y=genes_expre_p7_p8,
#                             by=c('gene_id'='gene_id'))
#
#expre_junto$tracking_id.x=NULL

#vroom_write(expre_junto, 
#           file = 'RNAseq_FPKM_1_to_8_merged.csv', 
#           delim = ',')

#hare un subset para q mi compu no muera

subset_FPKM <-  vroom(file = dir)

##Limpiar columna 'gene_id': quitar los puntos del gene_id (tomo parte del código de sol)

identificadores<- pull(subset_FPKM,gene_id)

identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])

subset_FPKM <- subset_FPKM %>% add_column(identificadores)


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

expre <-left_join(x = subset_FPKM,
                  y = myannot,
                  by = c('identificadores'='ensembl_gene_id'))

##Quedarme solo con los datos de tipo 'gene coding' (pues por ahora son los únicos que nos interesan)

protcod <- filter(expre, 
                       gene_biotype =='protein_coding')

valores_expre <- protcod %>% 
  dplyr::select(-gene_id,
                -gene_biotype,
                -identificadores,
                -percentage_gene_gc_content)

####  Discretizacion de los datos
mat_dis<-infotheo::discretize(t(valores_expre))

#Save matrix

vroom_write(mat_dis, 
        file = 'coexpression_matrix_AD.txt', 
       delim = ',')

###De aqui el siguiente sript es script_calculo_mi_parallel.R