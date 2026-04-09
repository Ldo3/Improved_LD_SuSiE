sim_2pop <- function(chrom, start, end, T_split, n_gwas, n_ref, J, h2, seed) {
  
  # Import necessary Python libraries
  stdpopsim <- reticulate::import("stdpopsim")
  msprime   <- reticulate::import("msprime")
  
  # Set R seed for downstream native R randomness
  set.seed(seed)
  
  # 1. Get Realistic Genetic Map from stdpopsim
  species <- stdpopsim$get_species("HomSap")
  contig  <- species$get_contig(chrom, left=as.integer(start), right=as.integer(end))
  
  # 2. Define Custom Demography
  demography <- msprime$Demography()
  demography$add_population(name="ANC", initial_size=20000L)
  demography$add_population(name="POP1", initial_size=100000L) # GWAS Cohort
  demography$add_population(name="POP2", initial_size=100000L) # LD Reference
  
  # reticulate translates R lists into Python lists seamlessly
  demography$add_population_split(time=as.integer(T_split), 
                                  derived=list("POP1", "POP2"), 
                                  ancestral="ANC")
  
  # 3. Run Simulation
  # Combine DTWF with standard coalescent
  models <- list(
    msprime$DiscreteTimeWrightFisher(duration=100L),
    msprime$StandardCoalescent()
  )
  
  ts <- msprime$sim_ancestry(
    samples = reticulate::dict(POP1 = as.integer(n_gwas), POP2 = as.integer(n_ref)),
    demography = demography,
    recombination_rate = contig$recombination_map,
    model = models,
    random_seed = as.integer(seed)
  )
  
  # Add mutations
  ts <- msprime$sim_mutations(ts, rate=1.29e-8, random_seed=1L)
  
  # 4. Extract Genotypes
  # Iterate to find numeric IDs for populations
  pop1_id <- NULL
  pop2_id <- NULL
  for (p in reticulate::iterate(ts$populations())) {
    if (!is.null(p$metadata) && p$metadata$name == "POP1") pop1_id <- p$id
    if (!is.null(p$metadata) && p$metadata$name == "POP2") pop2_id <- p$id
  }
  
  samples_gwas <- ts$samples(population=as.integer(pop1_id))
  samples_ref  <- ts$samples(population=as.integer(pop2_id))
  
  # Reticulate converts the numpy array automatically to an R matrix
  G_all <- ts$genotype_matrix() 
  
  # Helper to convert haplotypes to diploid dosage (0, 1, 2)
  to_diploid <- function(hap_idx) {
    # R is 1-indexed, Python hap_idx is 0-indexed
    idx <- as.integer(hap_idx) + 1 
    g_hap <- G_all[, idx, drop=FALSE]
    # Sum adjacent columns for diploid individuals
    g_dip <- g_hap[, seq(1, ncol(g_hap), by=2)] + g_hap[, seq(2, ncol(g_hap), by=2)]
    return(t(g_dip))
  }
  
  G_gwas <- to_diploid(samples_gwas) # Shape: N_gwas rows x N_snps cols
  G_ref  <- to_diploid(samples_ref)  # Shape: N_ref rows x N_snps cols
  
  # 5. Filter for Common Variants (MAF > 5% in both populations based on Python code)
  maf_gwas <- colMeans(G_gwas) / 2
  maf_ref  <- colMeans(G_ref) / 2
  
  common_mask <- (maf_gwas > 0.05) & (maf_gwas < 0.95) & (maf_ref > 0.05) & (maf_ref < 0.95)
  
  n_common <- sum(common_mask)
  if (n_common == 0) {
    stop("No common variants found. Try a larger region.")
  }
  
  G_gwas <- G_gwas[, common_mask, drop=FALSE]
  G_ref  <- G_ref[, common_mask, drop=FALSE]
  
  J_all <- max(ncol(G_gwas), J)
  
  # Random subset of J sequential SNPs
  max_start <- J_all - J
  start_idx <- if (max_start > 0) sample(0:max_start, 1) else 0
  
  # Adjusting for R's 1-based indexing for subsetting
  cols_to_keep <- (start_idx + 1):(start_idx + J)
  G  <- G_gwas[, cols_to_keep, drop=FALSE]
  G0 <- G_ref[, cols_to_keep, drop=FALSE]
  
  # 6. Generate summary stats
  # Center genotypes
  X <- sweep(G, 2, colMeans(G), "-")
  N <- nrow(X)
  
  G0_centered <- sweep(G0, 2, colMeans(G0), "-")
  
  # Calculate Covariance Matrices (rowvar=False in numpy = standard cov() in R)
  R0 <- cov(G0_centered)
  R0[is.na(R0)] <- 0 # Handle invariant SNPs
  
  R <- cov(X)
  R[is.na(R)] <- 0
  
  num_causal <- 2
  
  # Easy case
  first_causal_SNP <- sample(1:J, 1)
  second_causal_SNP <- which.min(R[first_causal_SNP, ]^2)
  causal_SNPs_loc <- c(first_causal_SNP, second_causal_SNP)
  
  beta_true <- rep(0, J)
  beta_true[causal_SNPs_loc] <- runif(num_causal) + 0.01
  
  # Scale heritability
  var_expl <- sum((X %*% beta_true)^2 / N)
  scale <- h2 / var_expl
  beta_true <- beta_true * sqrt(scale)
  
  # Simulate Phenotype y
  y <- (X %*% beta_true) + sqrt(1 - h2) * rnorm(N)
  
  # Note: Python's np.std() uses N denominator by default. R's sd() uses N-1. 
  # This exact equation replicates numpy's standard deviation scaling.
  y_std <- sqrt(sum((y - mean(y))^2) / N)
  y <- (y - mean(y)) / y_std
  
  # Calculate marginal effects
  v <- (t(X) %*% y) / N
  
  return(list(
    v = as.vector(v), 
    R = R, 
    R0 = R0, 
    beta_true = beta_true, 
    causal_SNPs_loc = causal_SNPs_loc
  ))
}