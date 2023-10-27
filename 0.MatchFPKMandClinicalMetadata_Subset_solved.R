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

#subset_counts <- counts[1:100,] #subset para darle chance a mi compu, en sefirot corro lo otro
subset_counts <- counts[1:100,-c(460,544,573)] # se eliminan gene ID, la copia del gene ID y dos muestras repetidas con _6 y _7 (se dejó el _0)

#nota:perdi 2 muestras (492_120515_6 y 492_120515_7)

#le quito los dos ultimos caracteres a los identificadores de geneID para matchearlos con 
#la metadata de RNASeq

colnames_no_num <- substr(colnames(subset_counts),1, nchar(colnames(subset_counts))-2)
colnames_no_num[1] <- "gene_id"  #como le quitamos los dos ultimos caracteres hay que recomponer el gene_id

# filtrar columnas específicas en la matriz subset_counts basándose en los valores
#de "specimenID" en RNA_seq_metadata y luego agregar una nueva columna "gene_id" en
#la primera posición del resultado.

subset_fpkm_matrix <- subset_counts[,(colnames_no_num %in% RNA_seq_metadata$specimenID)] %>% 
                  mutate(gene_id = subset_counts$gene_id, .before = 1)

RNA_seq_metadata <- RNA_seq_metadata[(RNA_seq_metadata$specimenID %in% colnames_no_num),c(19,1,17)] %>%  # 19 = specimenID 1 = individualID 17 = cogdx
                  arrange(match(specimenID, colnames_no_num[-1]))

#Recordando que
#1 NCI: No cognitive impairment (No impaired domains)
#2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
#3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
#4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)
#5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
#6 Other dementia: Other primary cause of dementia

#por lo cual no necesitamos el cogdx == 6

#libreria para cogdx == 1 (sanos)

cogdx1 <- RNA_seq_metadata %>% 
  filter(cogdx == 1)

#libreria para cogdx == c(2, 3) (MCI)

cogdx2_3 <- RNA_seq_metadata %>% 
  filter(cogdx == 2 | cogdx == 3)

#libreria para cogdx == c(4, 5) (AD)

cogdx4_5 <- RNA_seq_metadata %>% 
  filter(cogdx == 4 | cogdx == 5)

###########Filtrado final



#transversion 

t_counts <- subset_counts[,2:nrow(subset_counts)] %>% 
  t() %>%  #trans
  as.data.frame() %>% #lo vuelvo df
rownames_to_column(var = "specimenID") #hago que los nombres de filas sean una columna
