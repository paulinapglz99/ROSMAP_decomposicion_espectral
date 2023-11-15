#codigo para poner 0s en diagonal de una matriz de informacion mutua


#llamar a la red

matrix <- vroom(file = '/datos/rosmap/matriz_coexpre_allAD_11052023.txt') %>% 
  as.data.frame()

gene_names <- matrix %>% 
 pull(gene)

#debido a la estructura de nuestros datos, 
matrix <- matrix[,-1]

##Poner en 0 la diagonal ---

diag(matrix) <- 0

#volviendo a pegar la col ---

matrix <- data.frame(gene = gene_names, matrix)

#guardar ---

vroom_write(matrix, 
          file = '/datos/rosmap/matriz_coexpre_allAD_11052023_zero.txt', 
        delim = ',')
