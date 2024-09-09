#5.1 rds_to_matrix_parallel.R
#Script to convert RDS file to matrix and equal diagonal MI matrix to 0
#Script by Aidee Lashmi and paulinapglz.99@gmail.com

# Libraries --- ---
pacman::p_load("tidyverse", "Matrix", "purrr")

# Function to process a matrix --- ---
process_matrix <- function(input_path, output_path) {
  # Read the Mutual Information matrix
  mi_list <- readRDS(input_path)
  
  # Unnest the RDS
  coexpre_mat_unnested <- as.data.frame(do.call(cbind, mi_list)) %>% as.matrix()
  
  # Verify if the matrix is square
  if (nrow(coexpre_mat_unnested) != ncol(coexpre_mat_unnested)) {
    stop("Matrix is not square!")
  }
  
  # Test diagonal values
  IM_test <- coexpre_mat_unnested[, 5] %>% as.data.frame() # Adjust the column index if needed
  max_val <- max(IM_test)
  if (coexpre_mat_unnested[5, 5] != max_val) {
    warning("Diagonal value is not the highest in the column!")
  }
  
  # Set the diagonal to 0
  diag(coexpre_mat_unnested) <- 0
  
  # Keep only the upper triangle
  coexpre_mat_unnested[lower.tri(coexpre_mat_unnested, diag = TRUE)] <- NA
  
  # Save the processed matrix
  saveRDS(coexpre_mat_unnested, file = output_path)
}

# Define paths --- ---
input_paths <- c("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds",
                 "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds")

output_paths <- c("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds",
                  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds")

# Process matrices --- ---
mapply(process_matrix, input_paths, output_paths)

#END