import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def load_and_merge_data():
    expression = pd.read_csv("sig_exp.csv")
    peaks = pd.read_csv("sig_peaks.csv")
    peaks_coord = pd.read_csv("beddata.csv")
    genes_coord = pd.read_csv("gene_coordinates.csv")

    sig_peaks = pd.merge(peaks, peaks_coord, how='left', on='PeakID')
    sig_genes = pd.merge(expression, genes_coord, how='left', on="gene")
    return sig_peaks, sig_genes


def find_pax5_data(sig_genes):
    pax5_gene = sig_genes[sig_genes['gene'] == 'Pax5']
    pax5_start, pax5_end = int(pax5_gene.iloc[0]['start']), int(pax5_gene.iloc[0]['end'])
    return pax5_start, pax5_end


def filter_peaks_around_tss(sig_peaks, pax5_tss, upstream_range=1000000, downstream_range=1000000):
    sig_peaks = sig_peaks[sig_peaks["chrom"] == "chr4"]

    # Calculate the midpoint of each peak
    sig_peaks['midpoint'] = sig_peaks[['start', 'end']].mean(axis=1).astype(int)

    # Filter peaks upstream and downstream of the TSS
    filtered_peaks = sig_peaks[
        ((sig_peaks['midpoint'] >= (pax5_tss - upstream_range)) & (sig_peaks['midpoint'] <= pax5_tss)) |
        ((sig_peaks['midpoint'] >= pax5_tss) & (sig_peaks['midpoint'] <= (pax5_tss + downstream_range)))
        ]

    # Save filtered peaks to a BED file
    filtered_peaks[['PeakID','chrom', 'start', 'end']].to_csv('filtered_peaks_around_pax5.bed', sep='\t', index=False,
                                                     header=False)

    return filtered_peaks


def save_peaks_to_csv(near_pax5_peaks):
    near_pax5_peaks.to_csv("diff_peaks_around_Pax5.tsv", columns=['PeakID', 'chrom', 'start', 'end'],sep='\t', header = True,index=False)


def plot_peaks_near_pax5(near_pax5_peaks, genes_coordinates):
    fig, ax = plt.subplots(figsize=(10, 3))
    norm = mcolors.Normalize(vmin=near_pax5_peaks['logFC'].min(), vmax=near_pax5_peaks['logFC'].max())
    cmap = plt.get_cmap('Blues')
    gene_row = genes_coordinates[genes_coordinates['gene'] == 'Pax5'].iloc[0]
    ax.plot([gene_row['start'], gene_row['end']], [0.5, 0.5], color='lightgreen', linewidth=2, label='Pax5', marker='|')

    for _, row in near_pax5_peaks.iterrows():
        color = cmap(norm(row['logFC']))
        ax.plot([row['start'], row['end']], [1, 1], color=color, marker='|', markersize=10)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical')
    cbar.set_label('Log Fold Change')

    ax.set_title('Peaks Proximate to Pax5 Gene with Log Fold Change')
    ax.set_xlabel('Genomic Position')
    ax.set_ylim(0, 1.2)  # Adjust y-axis limits to compact space
    ax.set_yticks([])
    ax.legend()

    ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.tight_layout()
    plt.show()


def main():
    sig_peaks, sig_genes = load_and_merge_data()
    pax5_start, pax5_end = find_pax5_data(sig_genes)
    near_pax5_peaks = filter_peaks_around_tss(sig_peaks, pax5_start, 1000000, 1000000 )
    save_peaks_to_csv(near_pax5_peaks)
    genes_coordinates = pd.read_csv("gene_coordinates.csv")
    plot_peaks_near_pax5(near_pax5_peaks, genes_coordinates)


if __name__ == "__main__":
    main()
