import pandas as pd
import matplotlib.pyplot as plt

# Load data
atac_seq_data = pd.read_csv('atacseq.csv')
gene_expression_data = pd.read_csv('mmc2.csv', header=0, index_col=0)
gene_expression_data.index.name = 'Gene'

bed_columns = ['PeakID', 'chrom', 'start', 'end']
bed_data = pd.read_csv('beddata.txt', sep='\t', header=None, names=bed_columns)

# Merge ATAC-seq and BED data
merged_data = pd.merge(bed_data, atac_seq_data, on='PeakID')

# Define Pax5 region with upstream/downstream extension
upstream_extend = 200000
downstream_extend = 200000
pax5_peaks = merged_data[(merged_data['chrom'] == 'chr4') &
                         (merged_data['start'] <= 44710694 + downstream_extend) &
                         (merged_data['end'] >= 44524748 - upstream_extend)]

# Define subgroups
celltypes = atac_seq_data.columns
pre_b_cells = [
    "B.T1.Sp", "B.T2.Sp", "B.T3.Sp"
]
pro_b_cells = [
    "proB.CLP.BM", "proB.FrA.BM", "proB.FrBC.BM"
]
mature_b_cells = [
    "B.Fem.Sp", "B.Fo.Sp", "B.FrE.BM", "B.GC.CB.Sp", "B.GC.CC.Sp",
    "B.MZ.Sp", "B.PB.Sp", "B.PC.BM", "B.PC.Sp", "B.Sp", "B.mem.Sp", "B1b.PC"
]
other_cells = list(set(celltypes) - set(pre_b_cells + pro_b_cells + mature_b_cells + ['PeakID']))

# Extract data for each subgroup
subgroups = {
    'Pre-B Cells': pax5_peaks[pre_b_cells],
    'Pro-B Cells': pax5_peaks[pro_b_cells],
    'Mature B Cells': pax5_peaks[mature_b_cells],
    'Other Cells': pax5_peaks[other_cells]
}




import seaborn as sns


ridge_data = []
for group_name, group_data in subgroups.items():
    for peak, heights in zip(pax5_peaks['start'], group_data.values):
        ridge_data.append({'Height': heights.mean(), 'Subgroup': group_name})

ridge_df = pd.DataFrame(ridge_data)

sns.kdeplot(data=ridge_df, x='Height', hue='Subgroup', fill=True, common_norm=False)
plt.title('Density of Peak Heights by Subgroup')
plt.xlabel('Peak Height')
plt.ylabel('Density')
plt.show()

import matplotlib.pyplot as plt


pax5_start = 44524748
pax5_end = 44710694

interval_data = {
    "Pre-B Cells": subgroups["Pre-B Cells"],
    "Pro-B Cells": subgroups["Pro-B Cells"],
    "Mature B Cells": subgroups["Mature B Cells"],
    "Other Cells": subgroups["Other Cells"]
}

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 10), sharey=True)
axes = axes.flatten()

for i, (group_name, group_data) in enumerate(interval_data.items()):
    for _, row in pax5_peaks.iterrows():
        start = row['start']
        end = row['end']
        height = group_data.mean(axis=1).loc[row.name]  

        
        axes[i].plot([start, end], [height, height], c='blue', alpha=0.7)

    # Highlight Pax5 gene region
    axes[i].axvspan(pax5_start, pax5_end, color='red', alpha=0.2, label="Pax5 Region")
    axes[i].set_title(group_name)
    axes[i].set_xlabel("Genomic Coordinate (chr4)")
    axes[i].set_ylabel("Peak Height")
    axes[i].legend()

plt.tight_layout()
plt.suptitle("Peak Accessibility Intervals Across Subgroups", y=1.02, fontsize=16)
plt.show()
