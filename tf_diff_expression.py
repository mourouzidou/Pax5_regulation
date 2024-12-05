import os
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd
import matplotlib.pyplot as plt

# Enable automatic conversion between R and pandas
pandas2ri.activate()

# Install DESeq2 if not already installed
r('''.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2", update = FALSE)
''')

# Import DESeq2
deseq2 = importr('DESeq2')


def run_deseq2_per_cell_type(expression_file, output_dir):
    """
    Run DESeq2 analysis for each cell type against all others.
    Outputs all results without filtering by significance.
    """
    r(f'''
    library(DESeq2)
    library(data.table)

    # Load expression data
    expression_data <- read.csv("{expression_file}", row.names = 1)

    # Ensure data is integers
    expression_data <- round(expression_data)

    # Initialize list to store all results
    all_results_list <- list()

    # Loop over all cell types
    for (cell_type in colnames(expression_data)) {{
        # Create metadata
        metadata <- data.frame(
            row.names = colnames(expression_data),
            cell_type = ifelse(colnames(expression_data) == cell_type, cell_type, "Other")
        )

        # Create DESeq2 dataset
        dds <- DESeqDataSetFromMatrix(
            countData = expression_data,
            colData = metadata,
            design = ~ cell_type
        )

        # Run DESeq2
        dds <- DESeq(dds)

        # Extract results for the current cell type
        results <- results(dds, contrast = c("cell_type", cell_type, "Other"))
        results$gene <- rownames(results)

        # Save all results for the current cell type
        write.csv(as.data.frame(results), file = file.path("{output_dir}", paste0(cell_type, "_vs_all.csv")))

        # Add all results to results list
        all_results_list[[cell_type]] <- as.data.frame(results)
    }}

    # Combine all results across all cell types
    if (length(all_results_list) > 0) {{
        combined_all_results <- rbindlist(all_results_list, idcol = "cell_type")
        write.csv(as.data.frame(combined_all_results), file = file.path("{output_dir}", "combined_all_results.csv"))
    }} else {{
        cat("No results found across all cell types.")
    }}
    ''')


# Specify input files and parameters
expression_file = "tf_expr.csv"  # Input gene expression file
output_dir = "cell_type_deseq2_results"  # Output directory for DESeq2 results

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Run DESeq2
run_deseq2_per_cell_type(expression_file, output_dir)

# Load combined results in Python
combined_results_file = os.path.join(output_dir, "combined_all_results.csv")
combined_results = pd.read_csv(combined_results_file)

# Visualization: Heatmap of all results across cell types
heatmap_data = combined_results.pivot_table(
    index='gene', columns='cell_type', values='log2FoldChange', fill_value=0
)

plt.figure(figsize=(12, 8))
plt.imshow(heatmap_data, aspect='auto', cmap='coolwarm', interpolation='nearest')
plt.colorbar(label="Log2 Fold Change")
plt.xticks(range(heatmap_data.columns.shape[0]), heatmap_data.columns, rotation=90)
plt.yticks(range(heatmap_data.index.shape[0]), heatmap_data.index, fontsize=6)
plt.title("Differential Expression Across All Cell Types")
plt.tight_layout()

# Save the heatmap and CSV
heatmap_output_file = os.path.join(output_dir, "heatmap_cell_type_expression.pdf")
plt.savefig(heatmap_output_file)
plt.show()

print(f"Heatmap saved to {heatmap_output_file}")


