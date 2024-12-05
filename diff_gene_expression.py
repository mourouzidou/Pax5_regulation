import pandas as pd
from rpy2.robjects import pandas2ri, r, globalenv
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError

# pandas-R df conversion
pandas2ri.activate()# Replace with your file path

expression_df = pd.read_csv(gene_exp, index_col= 0)
celltypes = expression_df.columns
# B-cell ,Non-B-cell
b_cell_keywords = ["B.", "proB"]
b_cell_columns = [col for col in expression_df.columns if any(keyword in col for keyword in b_cell_keywords)]
non_b_cell_columns = list(set(expression_df.columns) - set(b_cell_columns))

b_cell_data = expression_df[b_cell_columns]
non_b_cell_data = expression_df[non_b_cell_columns]
combined_data = pd.concat([b_cell_data, non_b_cell_data], axis=1)

condition = [1] * len(b_cell_columns) + [0] * len(non_b_cell_columns)

ro.globalenv['combined_data'] = pandas2ri.py2rpy(combined_data)
ro.globalenv['condition'] = ro.IntVector(condition)
# R code 
r_code = """
library(limma)

# Create design matrix
design <- model.matrix(~ condition)

# Perform voom normalization
voom_data <- voom(as.matrix(combined_data), design)

# Fit the linear model
fit <- lmFit(voom_data, design)
fit <- eBayes(fit)

# Get the top differentially accessible peaks
results <- topTable(fit, coef=2, adjust="fdr", number=Inf)

# Perform one-sided test for positive fold changes
results <- results[results$logFC > 0, ]  # Keep only positive fold changes
results$P.Value <- results$P.Value / 2   # Adjust p-values for one-sided test
results
"""

try:
    results_r = r(r_code)
    if results_r.nrow > 0:
        results = pandas2ri.rpy2py_dataframe(results_r)
        significant_peaks = results[results['adj.P.Val'] < 0.05]
        print("Significant peaks with higher accessibility in B-cells:")
        print(significant_peaks)
        significant_peaks.to_csv("sig_exp.csv", index=False)
    else:
        print("No significant results found.")
except RRuntimeError as e:
    print("Error in R code:", e)
