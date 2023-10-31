# load libraries ---- 

library(tidyverse)
library(infotheo)
library(furrr)

#set timer ---- 

tempus <- Sys.time()

# read data ---- 



path1<-'Matriz_coexpre.txt'
solo_protcod <- vroom::vroom(file = path1)

# aqui hacemos el subsetting de la matriz de expresion ---- 

# solo_protcod <- solo_protcod %>% select() %>% filter()

# discretize ---- 

my_index <- 1:nrow(solo_protcod)
names(my_index) <- solo_protcod$identificadores[my_index]

plan(multicore, workers = 40)




datos_discretizados <- 
  future_map(.x = my_index, .f = function(i){
  solo_protcod[i,] %>% 
    dplyr::select(-gene_id,-gene_biotype,-identificadores,-percentage_gene_gc_content) %>% 
    unlist() %>% 
    discretize()
})

# calculate MI <- 

mi_mi <- 
  map(my_index, .f = function(i){
  ii = datos_discretizados[i]
  map(my_index, .f = function(j){
    jj = datos_discretizados[j]
    if(j>i){
    
      
      mutinformation(ii,jj)
      
    }
    
  })
  
})

#write out ---- 

write_rds(x =mi_mi, file = "/datos/rosmap/mi_rosmap_list.rds")

print(Sys.time() - tempus())