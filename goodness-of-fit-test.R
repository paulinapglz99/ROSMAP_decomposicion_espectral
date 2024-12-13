# Script de prueba de bondad de ajuste
#paulinapglz.99@gmail.com

pacman::p_load("dplyr", 
               "igraph", 
               "ggplot2", 
               "MASS", 
               "fitdistrplus", 
               "VGAM", 
               "poweRlaw")

#IMPORTANT NOTE === ====

# Null hypothesis (H0): Your data fit the power-law distribution well.
# p-value: The probability of getting as good (or better) a fit as observed, 
# under the assumption that the data actually follow a power law.

# If p> 0.05, i.e. if p is greater than the typical significance threshold (α=0.05), you do not reject H0.
# This means that your data can be modeled reasonably well with the tested distribution

#Functions --- ---

#Calculate degree distribution and freq table
calculate_degree_distribution <- function(graph, label) {
  degree_data <- degree(graph)
  degree_df <- data.frame(gene = names(degree_data), degree = degree_data)
  degree_freq <- table(degree_df$degree) %>% as.data.frame()
  colnames(degree_freq) <- c("degree", "Freq")
  degree_freq$degree <- as.numeric(as.character(degree_freq$degree))
  degree_freq$Prob <- degree_freq$Freq / sum(degree_freq$Freq)
  degree_freq <- degree_freq[order(degree_freq$degree, decreasing = FALSE), ]
  degree_freq$CumulativeDegree <- cumsum(degree_freq$Freq)
  degree_freq$logCumulativeDegree <- log10(degree_freq$CumulativeDegree)
  degree_freq$dx <- label
  degree_freq$CDF <- degree_freq$CumulativeDegree / sum(degree_freq$Freq)
  return(degree_freq)
}

#Function to fit distributions and calculate goodness of fitness

fit_distributions <- function(degree_freq) {
  results <- list()
  
  #Power-Law
  # pl_fit <- displ$new(degree_freq$degree)
  # est_pl <- estimate_xmin(pl_fit)
  # pl_fit$setXmin(est_pl)
  # p_value_powl <- bootstrap_p(pl_fit, no_of_sims = 500, threads = 4)
  # 
  # results$PowerLaw <- list(
  #   fit = pl_fit,
  #   p_value = p_value_powl
  # )
  # 
  # print("The p-value of the bootstrap for Power-Law distribution is:", results$PowerLaw$p_value, "\n")
  # 
  
  #Poisson
  poisson_fit <- fitdist(degree_freq$degree, "pois")
  lambda <- poisson_fit$estimate #Estimated lambda parameter
  #Kolmogorov-Smirnov test
  ks_test_pois <- ks.test(degree_freq$degree,"ppois", lambda)
  
  results$Poisson <- list(
    fit = poisson_fit,
    gof = gofstat(poisson_fit), 
    p_value = ks_test_pois$p.value
  )
  print("The p-value of the KS test for Poisson distribution is:", results$Poisson$p_value, "\n")
  
  # #Normal or gaussian dis
  # gaussian_fit <- fitdist(degree_freq$degree, "norm")
  # mean_est <- gaussian_fit$estimate["mean"]
  # sd_est <- gaussian_fit$estimate["sd"]
  # #Kolmogorov-Smirnov test
  # ks_test_norm <- ks.test(degree_freq$degree, "pnorm", mean_est, sd_est)
  # 
  # results$Gaussian <- list(
  #   fit = gaussian_fit,
  #   gof = gofstat(gaussian_fit), 
  #   p_value = ks_test_norm$p.value
  # )
  # 
  # cat("The p-value of the KS test for Gaussian distribution is:", results$Gaussian$p_value, "\n")
  # 
  # #Pareto dis
  # pareto_fit <- vglm(degree_freq$degree ~ 1, paretoff, trace = FALSE)
  # #Extract estimated parameters
  # shape <- Coef(pareto_fit)["shape"]
  # scale <- Coef(pareto_fit)["scale"]
  # 
  # #KS test
  # ks_test_pareto <- ks.test(degree_freq$degree, "ppareto", scale = scale, shape = shape)
  # 
  # results$Pareto <- list(
  #   fit = pareto_fit,
  #   gof = NULL, #The specific function for Pareto's gof may vary
  #   p_value =  ks_test_pareto$p.value
  #   )
  # 
  # cat("The p-value of the KS test for Pareto distribution is:", results$Pareto$p_value, "\n")
  # 
  # #Exponential
  # exp_fit <- fitdist(degree_freq$degree, "exp")
  # #Estimated lambda parameter
  # lambda_exp <- 1 / exp_fit$estimate["rate"]
  # 
  # #Kolmogorov-Smirnov test
  # ks_test_exp <- ks.test(degree_freq$degree, "pexp", rate = lambda_exp)
  # results$Exponential <- list(
  #   fit = exp_fit,
  #   gof = gofstat(exp_fit), 
  #   p_value = ks_test_exp$p.value
  # )
  # 
  # cat("The p-value of the KS test for exponential distribution is", results$Exponential$p_value, "\n")
  # 
  # #Negative binomial
  # negbin_fit <- fitdist(degree_freq$degree, "nbinom")
  # # Parámetros estimados
  # size_nb <- negbin_fit$estimate["size"]
  # prob_nb <- negbin_fit$estimate["prob"]
  # 
  # # Prueba de Kolmogorov-Smirnov
  # ks_test_nb <- ks.test(degree_freq$degree, "pnbinom", size = size_nb, prob = prob_nb)
  # 
  # results$NegativeBinomial <- list(
  #   fit = negbin_fit,
  #   gof = gofstat(negbin_fit), 
  #   p_value = ks_test_nb$p.value
  # )
  # 
  # cat("The p-value of the KS test for Negative binomial distribution is", results$NegativeBinomial$p_value, "\n")

  #Zero-inflated negative binomial
  # zeroinf_fit <- vglm(degree_freq$degree ~ 1, zinegbin, trace = FALSE)
  # results$ZeroInflatedNegativeBinomial <- list(
  #   fit = zeroinf_fit,
  #   gof = NULL # La función específica para gof de esta distribución también puede variar
  # )
  # 
  return(results)
}


#Read graph files
graph_files <- list(
  AD = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml",
  noAD = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml"
)

#List of graphs
graphLists <- lapply(graph_files, function(file) read_graph(file, format = "graphml"))

#Calculate degree dis for graphs
degree_distributions <- mapply(
  FUN = calculate_degree_distribution,
  graph = graphLists,
  label = names(graphLists),
  SIMPLIFY = FALSE
)

#Calculate again degree distributions and fit distributions
fits <- lapply(names(graphLists), function(label) {
  graph <- graphLists[[label]]
  degree_data <- degree(graph)
  degree_df <- data.frame(degree = as.numeric(degree_data)) %>%
    group_by(degree) %>%
    summarise(Freq = n())
  fit_results <- fit_distributions(degree_df)
  list(degree_df = degree_df, fit_results = fit_results)
})

#Names to degree distributions
names(degree_distributions) <- names(graphLists)

#See fits


