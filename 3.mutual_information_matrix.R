#Script que hace mutual information para datos de expresion RNAseq de ROSMAP

pacman::p_load('future', 
             'tidyverse', 
             'infotheo',
             'furrr', 
             'vroom')

#Set directory ---------

setwd(dir = '/datos/home/paulinapg/redesROSMAP/')

#set timer ---- 

tempus <- Sys.time()

#read data ----------

protcod <- vroom(file = 'protcod_MCI.txt')

mat_dis <- vroom(file = 'mat_dis_MCI.txt')

# indexing data --------------

my_index <- pull(protcod, 'gene_id') #de otro dataset not here pero que tenia en mi ambiente, ese esta en recap.R

#Quiero un vector que tenga los conteos de cada elemento de otro vector

my_index_i <- seq_along(my_index) #indice numero de my_index

names(my_index_i) <- my_index  #etiquetar con indices, le doy nombre a los elementos del vector

#set multicore plan --------

plan(multicore, workers = 40)

#Calculate mutual information ------

MI_MI <- future_map(.x = my_index_i, .f = function(k){
kk = mat_dis[k]
 map(.x = my_index_i, .f = function(m){
  mm = mat_dis[m]
   mutinformation(kk,mm)
 })
})

#write matrix ----

#write_rds(x = MI_MI,
#         file = "ROSMAP_RNAseq_MutualInfo_allMCI_matrix.rds") #this one did not work

saveRDS(MI_MI, "ROSMAP_RNAseq_MutualInfo_allMCI_matrix.rds")

#print time ----

print(Sys.time() - tempus)
