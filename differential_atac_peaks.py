import pandas as pd
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro

pandas2ri.activate()
atacseq_file = "atacseq.csv"  
atacseq_data = pd.read_csv(atacseq_file)

# Define B-cell and Non-B-cell columns
b_cell_columns = [
    "B.Fem.Sp", "B.Fo.Sp", "B.FrE.BM", "B.GC.CB.Sp", "B.GC.CC.Sp",
    "B.MZ.Sp", "B.PB.Sp", "B.PC.BM", "B.PC.Sp", "B.Sp",
    "proB.CLP.BM", "proB.FrA.BM", "proB.FrBC.BM"
]
non_b_cell_columns = list(set(atacseq_data.columns) - set(["PeakID", *b_cell_columns]))

peak_ids = atacseq_data["PeakID"]
b_cell_data = atacseq_data[b_cell_columns]
non_b_cell_data = atacseq_data[non_b_cell_columns]


combined_data = pd.concat([b_cell_data, non_b_cell_data], axis=1)
condition = [1] * len(b_cell_columns) + [0] * len(non_b_cell_columns)

r_combined_data = pandas2ri.py2rpy(combined_data)
r_condition = ro.IntVector(condition)

# R code for limma
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
results
"""

# Execute R code
r_env = ro.Environment()
r_env['combined_data'] = r_combined_data
r_env['condition'] = r_condition

results_r = ro.r(r_code, envir=r_env)

# Convert R results back to Python
results = pandas2ri.rpy2py(results_r)
results = results.reset_index()

# Save results
results.to_csv("differential_peaks.csv", index=False)

# Print significant peaks
significant_peaks = results[results['adj.P.Val'] < 0.05]
print("Significant peaks with higher accessibility in B-cells:")
print(significant_peaks)
