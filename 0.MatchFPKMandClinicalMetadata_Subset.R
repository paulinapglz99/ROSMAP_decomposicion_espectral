#Match Clinical.csv file to RNA Seq files and subset 

pacman::p_load("dplyr", 
               'vroom')
#Read files 

counts <- vroom( file = 'RNAseq_FPKM_1_to_8_merged.csv' ) #ROSMAP_RNAseq_FPKM_gene
metadata <- vroom( file = 'cli_bio_metadata.csv') #metadata para todos los ensayos

#Simplificar metadata a unicamente los ensayos de RNAseq en cerebro, y quitar
#los dos ultimos caracteres de las filas para matchear con specimenID

RNA_seq_metadata <- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain')
  
#hago un subset 

subset_counts <- counts[1:100,] #subset para darle chance a mi compu, en sefirot corro lo otro

#transversion 

t_counts <- subset_counts[,2:nrow(subset_counts)] %>% 
  t() %>%  #trans
  as.data.frame() %>% #lo vuelvo df
rownames_to_column(var = "specimenID_counts") #hago que los nombres de filas sean una columna

#de aqui en adelante no funciona jeje

t_counts$specimenID <- substring(t_counts$specimenID, 1, nchar(t_counts$specimenID) - 2)


# Realiza un filtro basado en la coincidencia parcial

