#Bayesian PRS by Mo Dehe, June 2024


# Install and load required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("medical-genomics-group/bayesR")
library(bayesR)
library(data.table)

# Load GWAS summary statistics
gwas_data <- fread("path/to/your/gwas_summary_stats.txt")

# Prepare GWAS data for BayesR
gwas_data <- gwas_data[, .(SNP = rsid, A1 = effect_allele, A2 = other_allele, 
                           N = sample_size, Z = effect_size / standard_error)]

# Write prepared GWAS data to a file
fwrite(gwas_data, "prepared_gwas_data.txt", sep = "\t")

# Set up BayesR parameters
genofile <- "path/to/your/genotype_data"  # PLINK binary format (.bed/.bim/.fam)
summaryfile <- "prepared_gwas_data.txt"
outfile <- "bayesR_output"

# Run BayesR
bayesR(genofile = genofile,
       summaryfile = summaryfile,
       outfile = outfile,
       n = max(gwas_data$N),  # Sample size
       maf = 0.01,            # Minor allele frequency threshold
       thread = 4,            # Number of threads to use
       S = c(0, 0.01, 0.1, 1),# Variance components
       pi = c(0.95, 0.02, 0.02, 0.01), # Prior probabilities for each component
       chain_length = 10000,  # Total number of MCMC iterations
       burn_in = 2000,        # Number of burn-in iterations
       thin = 5)              # Thinning interval

# Load BayesR results
effects <- fread(paste0(outfile, ".betas"))

# Calculate PRS
prs <- calculatePRS(genofile, effects)

# Save PRS
fwrite(data.table(IID = prs$IID, PRS = prs$PRS), 
       "bayesian_prs_results.txt", sep="\t")

# If you have phenotype data, you can evaluate the PRS
pheno <- fread("path/to/your/phenotype_data.txt")
combined_data <- merge(pheno, prs, by="IID")

# For binary outcome
logistic_model <- glm(phenotype ~ PRS, data=combined_data, family="binomial")
summary(logistic_model)

# For continuous outcome
linear_model <- lm(phenotype ~ PRS, data=combined_data)
summary(linear_model)

# Plot PRS distribution
library(ggplot2)
ggplot(combined_data, aes(x=PRS, fill=as.factor(phenotype))) +
  geom_density(alpha=0.5) +
  theme_minimal() +
  labs(title="PRS Distribution by Phenotype", 
       x="Polygenic Risk Score", 
       y="Density")

# Plot PRS vs phenotype (for continuous phenotype)
ggplot(combined_data, aes(x=PRS, y=phenotype)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_minimal() +
  labs(title="PRS vs Phenotype", 
       x="Polygenic Risk Score", 
       y="Phenotype")