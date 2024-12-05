import numpy as np
from mizani.transforms import log2_trans
from scipy.stats import zscore

from pax5_motifs import celltypes

# # Directories
# background_dir = "bg_fg/outputs_bg/aitac_bg_{i}/filter_predictions.npy"
#
# # Load all bg subsets
# bg_predictions = []
# for i in range(1, 31):
#     file_path = background_dir.format(i=i)
#     bg_predictions.append(np.load(file_path))
#
# bg_predictions_combined = np.vstack(bg_predictions)  # Shape: (n_bg_peaks, n_filters, n_cell_types)
# # print("Combined bg shape:", bg_predictions_combined.shape)
# # Load fg data
# fg_predictions = np.load("bg_fg/outputs_bg/aitac_fg/filter_predictions.npy")  # Shape: (n_fg_peaks, n_filters, n_cell_types)

# # Mean predictions across peaks
# fg_mean = fg_predictions.mean(axis=0)  # Shape: (n_filters, n_cell_types)
# bg_mean = bg_predictions_combined.mean(axis=0)  # Shape: (n_filters, n_cell_types)
#
# # Compute log-fold change
# log_fold_change = np.log2(fg_mean / (bg_mean + 1e-6))  # Add small value to avoid division by zero
# print("Log-fold change shape:", log_fold_change.shape)  # (n_filters, n_cell_types)
#
# # Threshold for significance
# significant_filters = np.where(np.abs(log_fold_change) > 1.0)  # Example threshold
#
# # List significant filters
# for filt, cell_type in zip(*significant_filters):
#     print(f"Filter {filt}, Cell Type {cell_type}, Log-Fold Change: {log_fold_change[filt, cell_type]}")

# import matplotlib.pyplot as plt
#
# # Example: Visualize top filters
# top_filters = np.argsort(-np.abs(log_fold_change.mean(axis=1)))[:10]
# for filt in top_filters:
#     plt.figure()
#     plt.plot(fg_mean[filt], label="Foreground (Pax5)")
#     plt.plot(bg_mean[filt], label="Background")
#     plt.title(f"Filter {filt} Activation")
#     plt.xlabel("Cell Type")
#     plt.ylabel("Mean Activation")
#     plt.legend()
#     plt.show()

# Load influence data
fg_influence = np.load("bg_fg/outputs_bg/aitac_fg/influence_by_OCR.npy")
bg_influence_combined = np.vstack([np.load(f"bg_fg/outputs_bg/aitac_bg_{i}/influence_by_OCR.npy") for i in range(1, 31)])
#
# # Mean influence scores
# fg_influence_mean = fg_influence.mean(axis=0)  # (n_filters,)
# bg_influence_mean = bg_influence_combined.mean(axis=0)  # (n_filters,)
#
# # Log-fold change in influence
# influence_log_fc = np.log2(fg_influence_mean / (bg_influence_mean + 1e-6))
# fg_normalized = zscore(fg_predictions, axis=0)  # Normalize FG along peaks
# bg_normalized = [zscore(bg, axis=0) for bg in bg_predictions_combined]  # Normalize BG subsets
# fg_mean = np.mean(fg_predictions, axis=0)
# bg_mean = np.mean(bg_predictions_combined, axis=0)
#
# log_fold_change = np.log2((fg_mean + 1e-6) / (bg_mean + 1e-6))  # Add small value to avoid division by zero
# from scipy.stats import ttest_ind
# significant_filters = []
# for f in range(300):
#     fg_values = fg_predictions[:, f, :]
#     bg_values = bg_predictions_combined[:, f, :]
#     t_stat, p_value = ttest_ind(fg_values.flatten(), bg_values.flatten(), equal_var=False)
#     if p_value < 0.05:
#         significant_filters.append(f)
#
# for f in significant_filters:
#     fg_cell_means = fg_predictions[:, f, :].mean(axis=0)
#     print(fg_cell_means)
#     bg_cell_means = bg_predictions_combined[:, f, :].mean(axis=0)
#     print()
#     cell_type_differences = fg_cell_means - bg_cell_means
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# # Load predictions (example arrays for FG and BG combined)
# # Assume fg_predictions and bg_predictions_combined are already loaded and normalized
#
# # Generate example data for demonstration
# np.random.seed(42)
# n_filters = 300
# n_cell_types = 81
# fg_predictions = np.random.rand(85, n_filters, n_cell_types)
# bg_predictions_combined = np.random.rand(30, 85, n_filters, n_cell_types)
#
# # Calculate mean FG and BG activations across peaks
# fg_mean = np.mean(fg_predictions, axis=0)  # Shape: (n_filters, n_cell_types)
# bg_mean = np.mean(bg_predictions_combined.reshape(-1, n_filters, n_cell_types), axis=0)
#
# # Calculate log-fold change (LFC) between FG and BG
# log_fold_change = np.log2((fg_mean + 1e-6) / (bg_mean + 1e-6))
# # Visualize LFC as a heatmap
# plt.figure(figsize=(12, 8))
# heatmap = sns.heatmap(log_fold_change, cmap="coolwarm", center=0, xticklabels=10, yticklabels=10, cbar_kws={"label": "LFC"})
# plt.title("Log-Fold Change (LFC) of Filter Activations: FG vs BG")
# plt.xlabel("Cell Types")
# plt.ylabel("Motif Filters")
# plt.tight_layout()
# plt.show()
#
# # Bar plot of top 10 filters with the highest mean LFC
# mean_lfc_per_filter = log_fold_change.mean(axis=1)
# top_filters = np.argsort(mean_lfc_per_filter)[-10:]
#
# plt.figure(figsize=(10, 6))
# plt.bar(range(10), mean_lfc_per_filter[top_filters], color='teal')
# plt.xticks(range(10), labels=top_filters, rotation=45)
# plt.title("Top 10 Filters with Highest Mean LFC (FG vs BG)")
# plt.xlabel("Filter Index")
# plt.ylabel("Mean LFC")
# plt.tight_layout()
# plt.show()
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# fg_predictions_combined = fg_predictions
#
# lfc_fg_vs_bg = np.log2(fg_predictions_combined + 1) - np.log2(bg_predictions_combined + 1)
# print("Log-Fold Change (LFC) shape:", lfc_fg_vs_bg.shape)
# mean_lfc_per_filter = np.mean(lfc_fg_vs_bg, axis=1)  # Mean across cell types
# top_20_filters_indices = np.argsort(mean_lfc_per_filter)[-20:]  # Indices of top 20 filters
# print("Top 20 filters by mean LFC:", top_20_filters_indices)
# mean_lfc_per_filter = np.mean(lfc_fg_vs_bg, axis=1)  # Mean across cell types
# top_20_filters_indices = np.argsort(mean_lfc_per_filter)[-20:]  # Indices of top 20 filters
# print("Top 20 filters by mean LFC:", top_20_filters_indices)
#
# mean_lfc_per_filter = np.mean(lfc_fg_vs_bg, axis=1)  # Mean across cell types
# top_20_filters_indices = np.argsort(mean_lfc_per_filter)[-20:]  # Indices of top 20 filters
# print("Top 20 filters by mean LFC:", top_20_filters_indices)
# import matplotlib.pyplot as plt
# plt.figure(figsize=(10, 8))
# plt.imshow(lfc_fg_vs_bg[top_20_filters_indices, :], aspect='auto', cmap='coolwarm', interpolation='nearest')
# plt.colorbar(label="LFC")
# plt.title("Log-Fold Change (LFC) of Top 20 Filters Across Cell Types")
# plt.xlabel("Cell Types")
# plt.ylabel("Top 20 Motif Filters")
# plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load FG and BG data
fg_influence_by_ocr = np.load("bg_fg/nofilter_fg/influence_by_OCR.npy")  # Shape: (85, 300)
bg_influence_by_ocr = np.stack([
    np.load(f"bg_fg/nofilter_bg/nofilter_bg_{i}/influence_by_OCR.npy")
    for i in range(1, 31)], axis=0)  # Shape: (30, 85, 300)

# Average BG influences across all subsets (mean over the 0th axis)
bg_mean_influence = np.mean(bg_influence_by_ocr, axis=0)  # Shape: (85, 300)

# Compute the mean influence across peaks for FG and BG
fg_mean_per_filter = np.mean(fg_influence_by_ocr, axis=0)  # Shape: (300,)
bg_mean_per_filter = np.mean(bg_mean_influence, axis=0)  # Shape: (300,)

# Calculate Log-Fold Change (LFC) between FG and BG
lfc_fg_vs_bg = np.log2((fg_mean_per_filter + 1e-8) / (bg_mean_per_filter + 1e-8))

# Identify top filters with the highest LFC in FG
top_filters_indices = np.argsort(lfc_fg_vs_bg)[::-1]  # Sort in descending order
top_20_filters = top_filters_indices[:20]  # Top 20 filters

# Visualization
plt.figure(figsize=(12, 6))
plt.bar(range(len(lfc_fg_vs_bg)), lfc_fg_vs_bg, color="gray", label="All Filters")
plt.bar(top_20_filters, lfc_fg_vs_bg[top_20_filters], color="darkmagenta", label="Top 20 Filters")
plt.axhline(0, color="black", linestyle="--", linewidth=0.8)
plt.xlabel("Filter Index")
plt.ylabel("Log-Fold Change (LFC)")
plt.title("Log-Fold Change of Filter Influence (FG vs. BG)")
plt.legend()
plt.tight_layout()
plt.show()

# Save the LFC results
np.save("lfc_fg_vs_bg.npy", lfc_fg_vs_bg)

# Highlight Top Filters
top_filters_df = pd.DataFrame({
    "Filter Index": top_20_filters,
    "FG Influence": fg_mean_per_filter[top_20_filters],
    "BG Influence": bg_mean_per_filter[top_20_filters],
    "LFC": lfc_fg_vs_bg[top_20_filters]
})
# print(top_filters_df)
# top_filters_df.to_csv("top_20_filters_fg_vs_bg.csv", index=False)


#___________________________________
# CELL WISE INFL PER OCR
#--------------------------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load FG and BG cellwise influence data
fg_cellwise_influence = np.load("bg_fg/nofilter_fg/cellwise_influence_by_OCR.npy")  # Shape: (85, 300, 81)
bg_cellwise_influence = np.stack([
    np.load(f"bg_fg/nofilter_bg/nofilter_bg_{i}/cellwise_influence_by_OCR.npy")
    for i in range(1, 31)], axis=0)  # Shape: (30, 85, 300, 81)

# Average BG influence across subsets
bg_mean_cellwise_influence = np.mean(bg_cellwise_influence, axis=0)  # Shape: (85, 300, 81)

# Mean influence across peaks for FG and BG
fg_mean_per_filter = np.mean(fg_cellwise_influence, axis=(0, 2))  # Shape: (300,) (mean over peaks and cell types)
bg_mean_per_filter = np.mean(bg_mean_cellwise_influence, axis=(0, 2))  # Shape: (300,)

# Log-Fold Change (LFC) between FG and BG for filters
lfc_fg_vs_bg = np.log2((fg_mean_per_filter + 1e-8) / (bg_mean_per_filter + 1e-8))

# Identify Top 20 Filters
top_filters_indices = np.argsort(lfc_fg_vs_bg)[::-1][:20]  # Indices of the top 20 filters

# Extract influence for Top 20 Filters across all cell types
fg_top_filters_cellwise = fg_cellwise_influence[:, top_filters_indices, :]  # Shape: (85, 20, 81)
bg_top_filters_cellwise = bg_mean_cellwise_influence[:, top_filters_indices, :]  # Shape: (85, 20, 81)

# Mean influence per cell type for FG and BG
fg_celltype_mean = np.mean(fg_top_filters_cellwise, axis=(0, 1))  # Shape: (81,)
bg_celltype_mean = np.mean(bg_top_filters_cellwise, axis=(0, 1))  # Shape: (81,)

# Load cell type names from atac.csv
atac = pd.read_csv("atacseq.csv")
cell_type_names = list(atac.columns)[1:]  # Assuming the column headers are the cell type names

# --- Visualization ---
plt.figure(figsize=(12, 6))

x = range(len(cell_type_names))  # Use cell type names as x-axis labels

plt.bar(x, fg_celltype_mean, width=0.4, label="FG (Foreground)", color="darkcyan", align="center")
plt.bar(x, bg_celltype_mean, width=0.4, label="BG (Background)", color="wheat", align="edge")
plt.xticks(x, cell_type_names, rotation=90)  # Replace x-axis indices with cell type names
plt.xlabel("Cell Types")
plt.ylabel("Mean Influence")
plt.title("Mean Influence by Cell Type for Top 20 Filters")
plt.legend()
plt.tight_layout()
# --- Save Plot as PNG ---
output_image_path = "top_20_filters_cellwise_influence_fg_vs_bg.png"
plt.savefig(output_image_path, dpi=300, bbox_inches='tight', format='png')  # Save as PNG
output_image_path_jpeg = "top_20_filters_cellwise_influence_fg_vs_bg.jpg"
plt.savefig(output_image_path_jpeg, dpi=300, bbox_inches='tight', format='jpeg')  # Save as JPEG



plt.show()

# Save Results
results = {
    "Top Filters": top_filters_indices,
    "FG Mean Influence by Cell Type": fg_celltype_mean,
    "BG Mean Influence by Cell Type": bg_celltype_mean,
    "LFC": lfc_fg_vs_bg[top_filters_indices],
}

np.save("top_20_filters_cellwise_influence_fg_vs_bg.npy", results)
