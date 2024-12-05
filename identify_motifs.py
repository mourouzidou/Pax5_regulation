import pandas as pd
import glob

all_motif_entries = []

for filename in glob.glob('*_fasta.txt'):
    file_id = filename.split('_')[0]

    # Process each file
    with open(filename, "r") as file:
        for line in file:
            if line.startswith('>'):
                parts = line.strip().split('_')
                chrom_info, offset_info = parts[0][1:], parts[2]  # Extract chromosome info and offset
                chrom, coords = chrom_info.split(':')
                start, end = map(int, coords.split('-'))
                offset = int(offset_info.split('=')[1].split()[0])  # Adjust if 'RC' present
                strand = 'RC' if 'RC' in offset_info else '+'

                # Adjust the start and end based on the offset and strand
                adjusted_start = start + offset if strand == '+' else start
                adjusted_end = end if strand == '+' else end - offset

                # Read the next line for the motif sequence
                motif_seq = next(file).strip()

                # Append a dictionary for each motif to the list
                all_motif_entries.append({
                    "Chromosome": chrom,
                    "Start": adjusted_start,
                    "End": adjusted_end,
                    "Motif": motif_seq,
                    "Strand": strand,
                    "SourceFile": file_id
                })

all_motif_data = pd.DataFrame(all_motif_entries)
print(all_motif_data.head())

genes_data = pd.read_csv('gene_coordinates.csv')
# all_motif_data.to_csv("combined_motif_locations.csv", index=False)

def check_overlap(motif_chrom, motif_start, motif_end, gene_chrom, gene_start, gene_end, threshold=1000):
    # Check if the motif and gene are on the same chromosome
    if motif_chrom != gene_chrom:
        return False
    # Check for overlap or proximity within the threshold
    if (motif_start <= gene_end + threshold) and (motif_end >= gene_start - threshold):
        return True
    return False

motif_data = pd.read_csv("combined_motif_locations.csv")

# List to store mapping results
mapping_results = []

# Iterate over each motif
for index, motif in motif_data.iterrows():
    # Check against each gene
    for _, gene in genes_data.iterrows():
        if check_overlap(motif['Chromosome'], motif['Start'], motif['End'],
                         gene['chrom'], gene['start'], gene['end']):
            mapping_results.append({
                'MotifID': motif['SourceFile'],
                'MotifStart': motif['Start'],
                'MotifEnd': motif['End'],
                'GeneName': gene['gene'],
                'GeneStart': gene['start'],
                'GeneEnd': gene['end'],
                'GeneID': gene['clusterId']
            })

# Create a DataFrame from the results
mapped_motifs_to_genes = pd.DataFrame(mapping_results)

# Save the results to a CSV file
mapped_motifs_to_genes.to_csv('mapped_motifs_to_genes.csv', index=False)

print(mapped_motifs_to_genes.head())
