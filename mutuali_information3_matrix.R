#Script que hace mutual information para datos de expresion RNAseq de ROSMAP

#set timer ---- 

tempus <- Sys.time()

#read data ----------

protcod <- vroom(file = 'protcod_AD.txt')

mat_dis <- vroom(file = 'mat_dis_AD.txt')

#hacer mis indices

my_index <- pull(protcod, 'gene_id') #de otro dataset not here pero que tenia en mi ambiente, ese esta en recap.R

#Quiero un vector que tenga los conteos de cada elemento de otro vector

my_index_i <- seq_along(my_index) #indice numero de my_index

names(my_index_i) <- my_index  #etiquetar con indices, le doy nombre a los elementos del vector

#antes de correr

plan(multicore, workers = 40)

#info mutua de matriz 

MI_AD <- map(.x = my_index_i, .f = function(k){
  
  kk = mat_dis[k]
  
  map(.x = my_index_i, .f = function(m){
    
    mm = mat_dis[m]
    
    mutinformation(kk,mm)
    
  })
})

#when done

write_rds(x =MI_AD, file = "/datos/rosmap/ROSMAP_RNAseq_MutualInfo_AD_matrix.rds")

#faroleando

print(Sys.time() - tempus())
