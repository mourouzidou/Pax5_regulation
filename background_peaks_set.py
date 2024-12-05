import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the ATAC-seq data
atacseq_file = "atacseq.csv"  # Replace with your file path
atacseq = pd.read_csv(atacseq_file)

# Remove the "PeakID" column for calculations
numeric_data = atacseq.iloc[:, 1:]  # Exclude PeakID for numerical operations

# Calculate median and standard deviation for each peak across all cell types
atacseq["median_accessibility"] = numeric_data.median(axis=1)
atacseq["std_accessibility"] = numeric_data.std(axis=1)

# Plot histograms of median and standard deviation
plt.figure(figsize=(12, 6))

# Median Accessibility Plot
plt.subplot(1, 2, 1)
plt.hist(atacseq["median_accessibility"], bins=30, color='skyblue', alpha=0.7)
plt.axvline(1.32, color='red', linestyle='dashed', label="Threshold: 1.32")
plt.title("Median Accessibility Distribution")
plt.xlabel("Median Accessibility")
plt.ylabel("Frequency")
plt.legend()

# Standard Deviation Plot
plt.subplot(1, 2, 2)
plt.hist(atacseq["std_accessibility"], bins=30, color='salmon', alpha=0.7)
plt.axvline(0.37, color='red', linestyle='dashed', label="Threshold: 0.37")
plt.title("Standard Deviation of Accessibility Distribution")
plt.xlabel("Standard Deviation")
plt.ylabel("Frequency")
plt.legend()

plt.tight_layout()
plt.show()

# Threshold justification:
# - The second "step" in the median accessibility plot is around 50,000, indicating peaks with low accessibility.
# - The second "step" in the standard deviation plot is around 0.37, indicating low variability.

# Filter for peaks that meet both threshold criteria
selected_peaks = atacseq[
    (atacseq["median_accessibility"] <= 1.32) &
    (atacseq["std_accessibility"] <= 0.37)
]
selected_peaks = selected_peaks.drop(columns = ["median_accessibility", "std_accessibility"])
print(selected_peaks.columns)
n_subsets = 30
random_subsets = []
for i in range(n_subsets):
    subset = selected_peaks.sample(n=85, replace=True, random_state=i)  # sampling with replacement
    random_subsets.append(subset)

print(random_subsets[0])

import numpy as np

# Generate random subsets with replacement and store them in a list
random_subsets = [
    selected_peaks.sample(n=85, replace=False, random_state=i).to_numpy()
    for i in range(n_subsets)
]

# Save all subsets as a single .npy file
output_file = "random_subsets_low_peaks.npy"
np.save(output_file, np.array(random_subsets))

print(f"Random subsets saved to {output_file}")


subsets = np.load("random_subsets_low_peaks.npy", allow_pickle=True)
print(subsets)
print(subsets.shape)


