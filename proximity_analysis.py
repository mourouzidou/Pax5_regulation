# import numpy as np
# import plotly.graph_objects as go
# import pandas as pd
# from matplotlib import pyplot as plt
#
# # Load your data
# data = pd.read_csv('mapped_motifs_to_genes.csv')
#
# # Create a figure
# fig = go.Figure()
#
# # Add gene traces
# genes = data.drop_duplicates(subset=['GeneName'])
# for _, g_row in genes.iterrows():
#     fig.add_trace(go.Scatter(
#         x=[g_row['GeneStart'], g_row['GeneEnd']],
#         y=[g_row['GeneName'], g_row['GeneName']],
#         mode='lines',
#         name=f"Gene: {g_row['GeneName']}",
#         text=f"{g_row['GeneName']} ({g_row['GeneStart']} - {g_row['GeneEnd']})",
#         line=dict(color='black', width=2),
#         hoverinfo='text'
#     ))
#
# # Add motif traces
# motif_ids = data['MotifID'].unique()
# colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(motif_ids)))
# color_dict = dict(zip(motif_ids, colors))
#
# for motif_id in motif_ids:
#     motif_data = data[data['MotifID'] == motif_id]
#     for _, m_row in motif_data.iterrows():
#         fig.add_trace(go.Scatter(
#             x=[m_row['MotifStart'], m_row['MotifEnd']],
#             y=[m_row['GeneName'], m_row['GeneName']],
#             mode='markers',
#             name=f"Motif from {motif_id}",
#             marker=dict(color=color_dict[motif_id], size=8),
#             text=f"Motif {motif_id} ({m_row['MotifStart']} - {m_row['MotifEnd']})",
#             hoverinfo='text'
#         ))
#
# # Update layout
# fig.update_layout(
#     title='Genomic Locations of Motifs and Genes by Source',
#     xaxis_title='Genomic Position',
#     yaxis=dict(title='Genes and Motifs', showgrid=False, zeroline=False),
#     showlegend=True
# )
#
# # Show figure
# fig.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data
significant_peaks = pd.read_csv("sig_peaks.csv", index_col=0)
bed_data = pd.read_csv("beddata.csv")  # File containing start positions of peaks

# Merge significant peaks with bed data to get genomic positions
merged_data = significant_peaks.merge(bed_data, on="PeakID")

# Define the range to filter peaks
pax5_start = 44524756
pax5_end = 44710487
buffer = 200000  # Buffer around PAX5 coordinates
filter_start = pax5_start - buffer
filter_end = pax5_end + buffer

# Filter for chromosome 4 within the specified range
chr4_data = merged_data[
    (merged_data["chrom"] == "chr4") &
    (merged_data["start"] >= filter_start) &
    (merged_data["start"] <= filter_end)
]

# Sort by start position within chromosome 4
chr4_data = chr4_data.sort_values(by="start")

# Load your ATAC-seq data
atacseq_data = pd.read_csv("atacseq.csv", index_col=0)

# Define columns for B-cell and Non-B-cell data
b_cell_columns = ["B.Fem.Sp", "B.Fo.Sp", "B.FrE.BM", "B.GC.CB.Sp", "B.GC.CC.Sp",
                  "B.MZ.Sp", "B.PB.Sp", "B.PC.BM", "B.PC.Sp", "B.Sp",
                  "proB.CLP.BM", "proB.FrA.BM", "proB.FrBC.BM"]
non_b_cell_columns = list(set(atacseq_data.columns) - set(b_cell_columns + ["PeakID"]))

# Extract and calculate the mean accessibility for each cell type
b_cell_means = atacseq_data[b_cell_columns].mean(axis=1)
non_b_cell_means = atacseq_data[non_b_cell_columns].mean(axis=1)

# Calculate log2 fold change
log2_fold_change = np.log2(b_cell_means + 1) - np.log2(non_b_cell_means + 1)

# Subset and reorder fold change data for chromosome 4
chr4_fold_change = log2_fold_change.loc[chr4_data['PeakID']]

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(
    chr4_fold_change.values.reshape(-1, 1),
    cmap="viridis",
    annot=False,
    cbar_kws={"label": "Log2 Fold Change (B-Cells vs. Non-B-Cells)"},
    yticklabels=False
)
positions = chr4_data["start"]

plt.axhspan(
    positions.searchsorted(pax5_start),
    positions.searchsorted(pax5_end),
    color="orange",
    alpha=0.16,
    label="PAX5 Region"
)
plt.title("Log2 Fold Change of Accessibility Near PAX5 on Chromosome 4", fontsize=14)
plt.xlabel("Log2 Fold Change", fontsize=12)
plt.ylabel("Peaks (Sorted by Genomic Position)", fontsize=12)
plt.tight_layout()

# Show the plot
plt.show()

