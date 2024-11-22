import pandas as pd
from rpy2.robjects import pandas2ri, r, globalenv
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError

# Enable the pandas-R dataframe conversion
pandas2ri.activate()

atacseq_file = "atacseq.csv"  
atacseq_data = pd.read_csv(atacseq_file)

# Define columns for B-cell and Non-B-cell data
b_cell_columns = ["B.Fem.Sp", "B.Fo.Sp", "B.FrE.BM", "B.GC.CB.Sp", "B.GC.CC.Sp",
                  "B.MZ.Sp", "B.PB.Sp", "B.PC.BM", "B.PC.Sp", "B.Sp",
                  "proB.CLP.BM", "proB.FrA.BM", "proB.FrBC.BM"]
non_b_cell_columns = list(set(atacseq_data.columns) - set(["PeakID", *b_cell_columns]))

# Extract relevant data
b_cell_data = atacseq_data[b_cell_columns]
non_b_cell_data = atacseq_data[non_b_cell_columns]

# Combine into a single matrix for analysis
combined_data = pd.concat([b_cell_data, non_b_cell_data], axis=1)
condition = [1] * len(b_cell_columns) + [0] * len(non_b_cell_columns)

# Convert to R data frames
globalenv['combined_data'] = pandas2ri.py2rpy(combined_data)
globalenv['condition'] = ro.IntVector(condition)

# R code to run in Python environment
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

# Try executing R code and handle potential errors
try:
    results_r = r(r_code)
    # Check if the R object is a dataframe and has rows
    if results_r.nrow > 0:
        results = pandas2ri.rpy2py_dataframe(results_r)
        # Filter results for significant peaks
        significant_peaks = results[results['adj.P.Val'] < 0.05]
        print("Significant peaks with higher accessibility in B-cells:")
        print(significant_peaks)
    else:
        print("No significant results found.")
except RRuntimeError as e:
    print("Error in R code:", e)
