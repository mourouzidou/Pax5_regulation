import pandas as pd
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro

# pandas-R dataframe conversion
pandas2ri.activate()

# Load gene expression data
gene_exp = "mmc2.csv"  # Replace with your file path

expression_df = pd.read_csv(gene_exp, index_col= 0)
celltypes = expression_df.columns
# Define B-cell and Non-B-cell columns
b_cell_keywords = ["B.", "proB"]
b_cell_columns = [col for col in expression_df.columns if any(keyword in col for keyword in b_cell_keywords)]
non_b_cell_columns = list(set(expression_df.columns) - set(b_cell_columns))

# Extract relevant data
b_cell_data = expression_df[b_cell_columns]
non_b_cell_data = expression_df[non_b_cell_columns]

# Combine into a single matrix
combined_data = pd.concat([b_cell_data, non_b_cell_data], axis=1)

# Create a condition vector (1 for B-cells, 0 for Non-B-cells)
condition = [1] * len(b_cell_columns) + [0] * len(non_b_cell_columns)

ro.globalenv['combined_data'] = pandas2ri.py2rpy(combined_data)
ro.globalenv['condition'] = ro.IntVector(condition)
# R code for limma
r_code = """
library(limma)

# Create design matrix
design <- model.matrix(~ condition)

# voom normalization
voom_data <- voom(as.matrix(combined_data), design)

# Fit the linear model
fit <- lmFit(voom_data, design)
fit <- eBayes(fit)

# Get the top differentially expressed genes
results <- topTable(fit, coef=2, adjust="fdr", number=Inf)
results
"""

# Execute R code and fetch results
results_r = ro.r(r_code)

# Convert R results back to Python
results = pandas2ri.rpy2py(results_r)
results = results.reset_index()

# Save results
results.to_csv("differential_exp.csv", index=False)

# Print significant diff expressed genes
significant_genes = results[results['adj.P.Val'] < 0.05]


