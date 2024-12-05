# import pandas as pd
# import numpy as np
#
# data = np.load("top_20_filters_cellwise_influence_fg_vs_bg.npy", allow_pickle=True).item()
# top_filters_indices = data['Top Filters']
# print(list(top_filters_indices))
#
# lfc = data['LFC']  # Log fold change values
# fg_mean_influence = data['FG Mean Influence by Cell Type']
# bg_mean_influence = data['BG Mean Influence by Cell Type']
#
# from Bio import motifs
#
# top_filters_indices = [217, 224, 257, 260, 89, 240, 112, 123, 275, 15, 94, 227, 165, 235, 167, 220, 173, 262, 78, 8]
#
#
# def extract_motifs(file_path, filter_indices):
#     with open(file_path) as f:
#         lines = f.readlines()
#
#     extracted_motifs = {}
#     for i, line in enumerate(lines):
#         if line.startswith("MOTIF"):
#             motif_id = int(line.split("filter")[1])
#             if motif_id in filter_indices:
#                 motif_data = lines[i:i + 21]
#                 extracted_motifs[motif_id] = "".join(motif_data)
#
#     return extracted_motifs
#
#
# # Extract FG motifs
# fg_motifs = extract_motifs("fg_motifs/filter_motifs_pwm.meme", top_filters_indices)
#
# # Extract BG motifs (combine results across subsets)
# bg_motifs = {}
# for i in range(1, 31):
#     subset_motifs = extract_motifs(f"bg_motifs/filter_motifs_pwm_bg_{i}.meme", top_filters_indices)
#     bg_motifs[f"bg_{i}"] = subset_motifs
#
# # Save extracted motifs
# with open("fg_top_motifs_pwm.meme", "w") as f:
#     for idx, motif in fg_motifs.items():
#         f.write(motif)
#
# with open("bg_top_motifs_pwm_combined.meme", "w") as f:
#     for subset, motifs in bg_motifs.items():
#         for idx, motif in motifs.items():
#             f.write(f">Subset: {subset}, Filter: {idx}\n{motif}")
# def convert_to_plain_text_pwm(input_file, output_file):
#     """
#     Extracts PWMs from a MEME file and converts them into plain text format for Tomtom.
#     Args:
#         input_file (str): Path to the MEME file.
#         output_file (str): Path to save the plain text PWM file.
#     """
#     with open(input_file, "r") as file:
#         lines = file.readlines()
#
#     with open(output_file, "w") as out:
#         motif_count = 0
#         for line in lines:
#             line = line.strip()
#             if line.startswith("MOTIF"):
#                 motif_count += 1
#                 out.write(f">Motif_{motif_count}\n")
#             elif line.startswith("letter-probability matrix"):
#                 continue  # Skip header lines
#             elif line:
#                 # Write PWM values directly
#                 out.write(f"{line}\n")
#
#     print(f"Converted to plain text PWM format and saved to {output_file}")
#
#
# # Convert the FG and BG files
# convert_to_plain_text_pwm("fg_top_motifs_pwm.meme", "fg_top_motifs_pwm_plain.txt")
# convert_to_plain_text_pwm("bg_top_motifs_pwm_combined.meme", "bg_top_motifs_pwm_combined_plain.txt")
input_file = "bg_top_motifs_pwm_combined_plain.txt"  # Replace with your input file name
output_file = "bg_motifs_pwm_cleaned.meme"  # Replace with your desired output file name

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    last_line_was_annotation = False  # Track if the last line was an annotation line
    for line in infile:
        # Check if the line starts with '>'
        if line.startswith(">"):
            if not last_line_was_annotation:  # Write only the first annotation in a row
                outfile.write(line)
            last_line_was_annotation = True
        else:
            outfile.write(line)  # Write non-annotation lines
            last_line_was_annotation = False

print(f"Cleaned motifs saved to {output_file}")
