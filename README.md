# Bayesian-PRS
Bayesian method for PRS calculation which is not very sensitive to base/target overlap


```markdown
# BayesR GWAS Analysis Pipeline

This repository contains an R script for running a BayesR analysis on GWAS summary statistics, calculating polygenic risk scores (PRS), and evaluating these scores with phenotype data.

## Requirements

Make sure you have the following packages installed:

- `devtools`
- `bayesR`
- `data.table`
- `ggplot2`

You can install these packages using the following commands:

```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("medical-genomics-group/bayesR")
install.packages(c("data.table", "ggplot2"))
```

## Steps

### 1. Load and Prepare GWAS Data

First, load the GWAS summary statistics and prepare the data for BayesR:

```R
library(bayesR)
library(data.table)

# Load GWAS summary statistics
gwas_data <- fread("path/to/your/gwas_summary_stats.txt")

# Prepare GWAS data for BayesR
gwas_data <- gwas_data[, .(SNP = rsid, A1 = effect_allele, A2 = other_allele, 
                           N = sample_size, Z = effect_size / standard_error)]

# Write prepared GWAS data to a file
fwrite(gwas_data, "prepared_gwas_data.txt", sep = "\t")
```

### 2. Set Up BayesR Parameters and Run Analysis

Specify the genotype file and run the BayesR analysis:

```R
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
```

### 3. Calculate Polygenic Risk Scores (PRS)

Load the BayesR results and calculate the PRS:

```R
# Load BayesR results
effects <- fread(paste0(outfile, ".betas"))

# Calculate PRS
prs <- calculatePRS(genofile, effects)

# Save PRS
fwrite(data.table(IID = prs$IID, PRS = prs$PRS), 
       "bayesian_prs_results.txt", sep="\t")
```

### 4. Evaluate PRS with Phenotype Data

If you have phenotype data, you can merge it with the PRS and perform logistic or linear regression analysis:

```R
# If you have phenotype data, you can evaluate the PRS
pheno <- fread("path/to/your/phenotype_data.txt")
combined_data <- merge(pheno, prs, by="IID")

# For binary outcome
logistic_model <- glm(phenotype ~ PRS, data=combined_data, family="binomial")
summary(logistic_model)

# For continuous outcome
linear_model <- lm(phenotype ~ PRS, data=combined_data)
summary(linear_model)
```

### 5. Visualize PRS Distribution and Relationship with Phenotype

Use `ggplot2` to visualize the distribution of PRS and its relationship with the phenotype:

```R
library(ggplot2)

# Plot PRS distribution
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
```
