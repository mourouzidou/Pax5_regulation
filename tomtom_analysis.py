import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load Tomtom FG and BG results
fg_tomtom = pd.read_csv("tomtom_fg.tsv", sep="\t", comment="#")
bg_tomtom = pd.read_csv("tomtom_bg.tsv", sep="\t", comment="#")

# Set significance threshold
threshold = 0.05
fg_tomtom = pd.read_csv("tomtom_fg.tsv", sep="\t", comment="#")

# # Set threshold for E-value
# evalue_threshold = 1.0
#
# # Filter motifs based on E-value
# fg_significant_evalue = fg_tomtom[fg_tomtom["E-value"] < evalue_threshold]
#
# # Sort motifs by E-value
# fg_significant_evalue = fg_significant_evalue.sort_values(by="E-value")
#
# # Plot -log10(E-value) for enriched FG motifs
# plt.figure(figsize=(12, 8))
# plt.barh(
#     fg_significant_evalue["Query_ID"],
#     -np.log10(fg_significant_evalue["E-value"]),
#     color="dodgerblue"
# )
# plt.xlabel("-log10(E-value)")
# plt.ylabel("Motif ID")
# plt.title("Enriched FG Motifs (E-value < 1.0)")
# plt.gca().invert_yaxis()  # Invert y-axis for better readability
# plt.tight_layout()
# plt.show()

import pandas as pd

# Load the Tomtom results file
fg_tomtom_hocomoco = pd.read_csv("tomtom_fg_hocomoco.tsv", sep="\t", comment="#")

significant_motifs = fg_tomtom_hocomoco[fg_tomtom_hocomoco["q-value"] < 0.05]

target_ids = list(significant_motifs["Target_ID"].unique())
# print("Unique Target_IDs:", len(target_ids))
#
# a_quality_tfbs = [id for id in target_ids if id.split('.')[-1] =='A']
# print(a_quality_tfbs)
print(target_ids)
expressed_tfs = ['SPIB', 'IRF4', 'IRF5', 'IRF8', 'LYL']
tfs = list(map(lambda x : x.split(".")[0], target_ids))

filtered_motifs = [id for id in target_ids if any(id.startswith(gene) for gene in expressed_tfs)]
print(filtered_motifs)



#
# # fg_tomtom_hocomoco = fg_tomtom_hocomoco[fg_tomtom_hocomoco["Target_ID"].isin(a_quality_tfbs)]
# # Filter for significant motifs (optional if already filtered)
# threshold = 0.05  # Adjust if needed
# fg_tomtom_hocomoco = fg_tomtom_hocomoco[fg_tomtom_hocomoco["p-value"] < threshold]
# #
# # Identify the best match for each Query_ID based on the lowest p-value
# best_matches = fg_tomtom_hocomoco.loc[
#     fg_tomtom_hocomoco.groupby("Query_ID")["p-value"].idxmin()
# ]
#
# # # # Save the best matches to a new file
# # # best_matches.to_csv("fg_best_matches_hocomoco.tsv", sep="\t", index=False)
# #
# # Print the top results
# print(best_matches.head())
#
# # Extract Target_IDs for further analysis
# best_target_ids = best_matches["Target_ID"].unique()
# print("Best Target IDs:", best_target_ids)
#
# best_tfs = list(map(lambda x : x.split(".")[0], best_target_ids))
# print(set(best_tfs))
#
#
# expr = pd.read_csv("mmc2.csv")
# expr['GeneID'] =  expr["GeneID"].apply(lambda x: str(x).upper())
# expr = expr[expr['GeneID'].isin(tfs)]
# expr.to_csv('tf_expr.csv', index=False)
#
# b_cells = ['B.Fem.Sp', 'B.Fo.Sp', 'B.FrE.BM', 'B.GC.CB.Sp', 'B.GC.CC.Sp', 'B.MZ.Sp', 'B.PB.Sp', 'B.PC.BM',
#                       'B.PC.Sp', 'B.Sp',
#                       'B.T1.Sp', 'B.T2.Sp', 'B.T3.Sp', 'B.mem.Sp', 'B1b.PC']

#
# import os
# from rpy2.robjects import r, pandas2ri
# from rpy2.robjects.packages import importr
# import pandas as pd
#
# # Enable automatic conversion between R and pandas
# pandas2ri.activate()
#
# # Install DESeq2 if not already installed
# r('''
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2", update = FALSE)
# ''')
#
# # Import DESeq2
# deseq2 = importr('DESeq2')
#
# # Function to run DESeq2 and filter TFs
# def run_deseq2(expression_file, tf_list, output_file):
#     """
#     Run DESeq2 analysis and filter for transcription factors.
#     """
#     r(f'''
#     library(DESeq2)
#
#     # Load expression data
#     expression_data <- read.csv("{expression_file}", row.names = 1)
#
#     # Create metadata
#     metadata <- data.frame(
#         row.names = colnames(expression_data),
#         cell_type = c(rep("B_cell", 3), rep("Non_B_cell", ncol(expression_data) - 3))  # Adjust as needed
#     )
#
#     # Create DESeq2 dataset
#     dds <- DESeqDataSetFromMatrix(
#         countData = expression_data,
#         colData = metadata,
#         design = ~ cell_type
#     )
#
#     # Run DESeq2
#     dds <- DESeq(dds)
#
#     # Extract results for B cells vs Non-B cells
#     results <- results(dds, contrast = c("cell_type", "B_cell", "Non_B_cell"))
#     results$gene <- rownames(results)
#
#     # Filter significant genes (adjusted p-value < 0.05)
#     significant_results <- results[results$padj < 0.05, ]
#     significant_results <- significant_results[order(significant_results$log2FoldChange, decreasing = TRUE), ]
#
#     # Save full results
#     write.csv(as.data.frame(significant_results), file = "deseq2_results.csv")
#
#     # Filter for transcription factors
#     tfs <- c({", ".join(f'"{tf}"' for tf in tf_list)})
#     tf_results <- significant_results[significant_results$gene %in% tfs, ]
#     write.csv(as.data.frame(tf_results), file = "{output_file}")
#     ''')
#
# # Specify input files and parameters
# expression_file = "mmc2.csv"  # Input gene expression file
# tf_list = tfs  # Replace with your TFs of interest
# output_file = "tf_deseq2_results.csv"  # Output file for filtered TF results
#
# # Run DESeq2
# run_deseq2(expression_file, tf_list, output_file)
#
# # Load and visualize the results in Python
# tf_results = pd.read_csv(output_file)
# print("Filtered TF Results:")
# print(tf_results)
#
# # Visualization
# import matplotlib.pyplot as plt
#
# plt.figure(figsize=(10, 6))
# plt.barh(tf_results['gene'], tf_results['log2FoldChange'], color='steelblue')
# plt.xlabel("Log2 Fold Change")
# plt.ylabel("Transcription Factors")
# plt.title("Differential Expression of TFs in B Cells")
# plt.tight_layout()
# plt.show()
