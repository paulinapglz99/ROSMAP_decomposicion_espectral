#Script que hace mutual information para datos de expresion RNAseq de ROSMAP

pacman::p_load('future', 
               'tidyverse', 
               'infotheo',
               'furrr', 
               'vroom'
               )
#set directory ---

setwd(dir = '/datos/home/paulinapg/redesROSMAP/')

#set timer ---- 

tempus <- Sys.time()

#read data ----------

protcod <- vroom(file = 'protcod_AD.txt')  #this file contains 

mat_dis <- vroom(file = 'mat_dis_AD.txt')  #this file contains a coexpression discretized matrix 

#indexing for paralell ----

my_index <- pull(protcod, 'gene_id') 

#Make a vector with counts for every element 

my_index_i <- seq_along(my_index) #index number for my_index (indice del indice)

names(my_index_i) <- my_index  #etiquetar con indices, le doy nombre a los elementos del vector

#set multicore plan for paralell -------

plan(multicore, workers = 40)

#calculate mutual information ------

MI_MI <- future_map(.x = my_index_i, .f = function(k){
  
  kk = mat_dis[k]
  
  map(.x = my_index_i, .f = function(m){
    
    mm = mat_dis[m]
    
    mutinformation(kk,mm)
    
  })
})

#when done, write matrix in rds format ------

saveRDS(MI_MI, "ROSMAP_RNAseq_MutualInfo_allAD_matrix.rds")

#faroleanding ----

print(Sys.time() - tempus())
