import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

# Load the results
output_dir = "cell_type_deseq2_results"
combined_results_file = os.path.join(output_dir, "combined_all_results.csv")
combined_results = pd.read_csv(combined_results_file)

# Pivot the data for heatmap
heatmap_data = combined_results.pivot_table(
    index='gene', columns='cell_type', values='log2FoldChange', fill_value=0
)

# Emphasize upregulation (set negative values to 0)
heatmap_data = heatmap_data.clip(lower=0)

# Hierarchical clustering for genes only
row_linkage = linkage(heatmap_data, method="weighted")
row_order = leaves_list(row_linkage)
#
# # Reorder rows based on clustering, keep cell types in original order
# heatmap_data = heatmap_data.iloc[row_order, :]
#
# # Plot the heatmap
# plt.figure(figsize=(14, 10))
# sns.heatmap(
#     heatmap_data,
#     cmap="Reds",  # Use a red-only colormap to emphasize upregulation
#     cbar_kws={"label": "Log2 Fold Change (Upregulation)"},
#     linewidths=0.5
# )
# plt.title("Differential Expression Clustered by Genes (Upregulation Focus)")
# plt.xlabel("Cell Types")
# plt.ylabel("Genes")
# plt.tight_layout()
#
# # Save the plot and clustered data
# clustered_heatmap_file = os.path.join(output_dir, "clustered_heatmap_genes_only.pdf")
# plt.savefig(clustered_heatmap_file)
# plt.show()
#
# # Save the clustered data
# clustered_data_file = os.path.join(output_dir, "clustered_heatmap_genes_only.csv")
# heatmap_data.to_csv(clustered_data_file)
#
# print(f"Clustered heatmap saved to {clustered_heatmap_file}")
# print(f"Clustered data saved to {clustered_data_file}")


start = 44524748 -1_000_000
end = 44710694 + 1_000_000
print(start)
print(end-start)


