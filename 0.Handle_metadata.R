#libraries  ----- 

pacman::p_load("dplyr", 
               'vroom')

#Handle metadata ---------------

clinical_metadata <- vroom(file = 'ROSMAP_clinical.csv')

biospecimen_metadata <- vroom(file = 'ROSMAP_biospecimen_metadata.csv')

#merge metadata in one

cli_bio_metadata <- left_join(x = clinical_metadata, 
                              y = biospecimen_metadata, 
                              by = "individualID")

#save table for later

#vroom_write(cli_bio_metadata,
#           file = "cli_bio_metadata.csv", 
#           delim = ",")

#read expression data

FPKM_p1_p6 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv')
FPKM_p7_p8 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv')

#merge expression data

FPKM_p1_p8<-left_join(x=FPKM_p1_p6,
                      y=FPKM_p7_p8,
                      by=c('gene_id'='gene_id'))

#Match Clinical.csv file to RNA Seq files and subset -----------------------

#Here we handle with inopportune data like the copy of the gene ID col
#and two repeated samples with _6 and _7 are deleted (_0 was left)

#merge both plates

FPKM_p1_p8 <- FPKM_p1_p8 %>% 
  dplyr::select(-c("tracking_id.x", "tracking_id.y", #remove inopportune data 
            "492_120515_6", "492_120515_7")) 

#note: I lost 2 samples (492_120515_6 and 492_120515_7)

#save table for later

#vroom_write(FPKM_p1_p8,
#           file = "RNAseq_FPKM_1_to_8_merged.csv", 
#           delim = ",")
