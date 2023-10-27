# radian

#Match Clinical.csv file to RNA Seq files and subset 

pacman::p_load("dplyr", 
               'vroom', 
               'tibble')
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

colnames_no_num <- substr(colnames(subset_counts),1, nchar(colnames(subset_counts))-2)
colnames_no_num[1] <- "gene_id"

subset_fpkm_matrix <- subset_counts[,(colnames_no_num %in% RNA_seq_metadata$specimenID)] %>% 
                  mutate(gene_id = subset_counts$gene_id, .before = 1)

RNA_seq_metadata <- RNA_seq_metadata[(RNA_seq_metadata$specimenID %in% colnames_no_num),c(19,1,17)] %>%  # 19 = specimenID 1 = individualID 17 = cogdx
                  arrange(match(specimenID, colnames_no_num[-1]))


RNA_seq_metadata$cogdx == 4

el4 <- subset_fpkm_matrix[,c(TRUE,RNA_seq_metadata$cogdx == 4)]

colnames(el4)

#transversion 

t_counts <- subset_counts[,2:nrow(subset_counts)] %>% 
  t() %>%  #trans
  as.data.frame() %>% #lo vuelvo df
rownames_to_column(var = "specimenID") #hago que los nombres de filas sean una columna

#le quito los ultimos dos caracteres a la columna

t_counts$specimenID <- substr(t_counts$specimenID,1,
                                     nchar(t_counts$specimenID)-2)

# Realiza un filtro basado en la coincidencia de specimenID

#esto lo quito despues

specimenID_counts <- pull(t_counts, 
                   specimenID)

specimenID_metadata <- pull(RNA_seq_metadata, 
                           specimenID) 

coincidencias <- intersect(specimenID_counts,   #este es un vector que tiene todos los specimenIDs que tienen una anotacion en la metadata
                           specimenID_metadata)

####

joined <- left_join(t_counts,
                    RNA_seq_metadata, 
                    by = coincidencias)
              

