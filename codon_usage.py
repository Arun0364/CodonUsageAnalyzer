# ======================================================
# 🧬 Codon Usage Analyzer
# Author: Arun
#
# Description:
#   - Reads one or two DNA sequences from FASTA files
#   - Counts codons and calculates codon usage frequencies (%)
#   - Exports results as CSV
#   - Visualizes codon usage as:
#       1. Bar plots (per sequence)
#       2. Codon heatmaps (per sequence)
#       3. Side-by-side comparison plots (two sequences)
#       4. Difference heatmaps (Species1 - Species2)
#
# Usage:
#   python codon_usage.py
#   -> Script will prompt for FASTA file paths and labels
#
# Outputs:
#   Results are saved in `outputs/` folder:
#       - codon_usage_<species>.csv
#       - codon_plot_<species>.png
#       - codon_heatmap_<species>.png
#       - codon_comparison.png (if two sequences provided)
#       - codon_diff_heatmap.png (if two sequences provided)
# ======================================================

from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ---------------------------
# Sequence Processing
# ---------------------------

def read_fasta(file_path):
    """
    Read the first sequence from a FASTA file.
    
    Parameters:
        file_path (str): Path to the FASTA file
    
    Returns:
        str: DNA sequence in uppercase letters (A, T, G, C)
    """
    for record in SeqIO.parse(file_path, "fasta"):
        return str(record.seq).upper()
    return ""

def count_codons(sequence):
    """
    Count codons (non-overlapping triplets) in a DNA sequence.
    
    Parameters:
        sequence (str): DNA sequence
    
    Returns:
        dict: Dictionary of codon -> count
    """
    codon_counts = defaultdict(int)
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            codon_counts[codon] += 1
    return codon_counts

def calculate_frequencies(codon_counts):
    """
    Convert raw codon counts into relative frequencies (%).
    
    Parameters:
        codon_counts (dict): Dictionary of codon -> count
    
    Returns:
        dict: Dictionary of codon -> frequency (% of total codons)
    """
    total_codons = sum(codon_counts.values())
    if total_codons == 0:
        return {codon: 0.0 for codon in codon_counts}
    return {codon: (count / total_codons) * 100 for codon, count in codon_counts.items()}

# ---------------------------
# Visualization
# ---------------------------

def plot_codon_barplot(counts, freqs, label, output_path):
    """
    Generate a barplot of codon usage frequencies for one species.
    """
    df = pd.DataFrame({
        "Codon": list(counts.keys()),
        "Count": list(counts.values()),
        "Frequency": [freqs[c] for c in counts.keys()]
    }).sort_values(by="Codon")

    plt.figure(figsize=(14,6))
    sns.barplot(x="Codon", y="Frequency", hue="Codon", data=df,
                palette="viridis", legend=False)
    plt.xticks(rotation=90)
    plt.title(f"Codon Usage Frequency (%) - {label}")
    plt.ylabel("Frequency (%)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_codon_heatmap(freqs, label, output_path):
    """
    Generate a codon usage heatmap arranged in codon table order.
    """
    bases = ["U", "C", "A", "G"]
    codon_grid = []
    for b1 in bases:
        row = []
        for b2 in bases:
            for b3 in bases:
                codon = (b1+b2+b3).replace("U","T")
                row.append(freqs.get(codon, 0.0))
        codon_grid.append(row)

    col_labels = [b2+b3 for b2 in bases for b3 in bases]
    df = pd.DataFrame(codon_grid, index=bases, columns=col_labels)

    plt.figure(figsize=(16,6))
    sns.heatmap(df, annot=True, fmt=".1f", cmap="YlOrRd", 
                cbar_kws={'label': 'Frequency (%)'})
    plt.title(f"Codon Usage Heatmap (%) - {label}")
    plt.xlabel("Second+Third Base")
    plt.ylabel("First Base")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def compare_codon_usage(freqs1, freqs2, label1, label2, output_path):
    """
    Compare codon usage between two sequences using side-by-side barplots.
    """
    codons = sorted(set(freqs1.keys()) | set(freqs2.keys()))
    df = pd.DataFrame({
        "Codon": codons,
        label1: [freqs1.get(c, 0) for c in codons],
        label2: [freqs2.get(c, 0) for c in codons]
    })
    df_melted = df.melt(id_vars="Codon", var_name="Species", value_name="Frequency")

    plt.figure(figsize=(16,7))
    sns.barplot(x="Codon", y="Frequency", hue="Species", data=df_melted, palette="viridis")
    plt.xticks(rotation=90)
    plt.title(f"Codon Usage Comparison: {label1} vs {label2}")
    plt.ylabel("Frequency (%)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def codon_difference_heatmap(freqs1, freqs2, label1, label2, output_path):
    """
    Generate a heatmap showing codon usage differences (Species1 - Species2).
    """
    bases = ["U", "C", "A", "G"]
    codon_grid = []
    for b1 in bases:
        row = []
        for b2 in bases:
            for b3 in bases:
                codon = (b1+b2+b3).replace("U","T")
                diff = freqs1.get(codon, 0) - freqs2.get(codon, 0)
                row.append(diff)
        codon_grid.append(row)

    col_labels = [b2+b3 for b2 in bases for b3 in bases]
    df = pd.DataFrame(codon_grid, index=bases, columns=col_labels)

    plt.figure(figsize=(16,6))
    sns.heatmap(df, annot=True, fmt=".1f", cmap="coolwarm", center=0,
                cbar_kws={'label': f'Difference (%) {label1}-{label2}'})
    plt.title(f"Codon Usage Difference Heatmap: {label1} vs {label2}")
    plt.xlabel("Second+Third Base")
    plt.ylabel("First Base")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# ---------------------------
# Export
# ---------------------------

def export_to_csv(counts, freqs, label, output_path):
    """
    Export codon counts and frequencies to a CSV file.
    """
    df = pd.DataFrame({
        "Codon": list(counts.keys()),
        "Count": list(counts.values()),
        "Frequency (%)": [freqs[c] for c in counts.keys()]
    }).sort_values(by="Codon")
    df.to_csv(output_path, index=False)

# ---------------------------
# Main Program
# ---------------------------

if __name__ == "__main__":
    os.makedirs("outputs", exist_ok=True)

    # === User Input ===
    fasta1 = input("Enter path to first FASTA file (e.g., sequences/species1.fasta): ").strip()
    label1 = input("Enter label for first sequence (e.g., Human): ").strip()
    
    fasta2 = input("Enter path to second FASTA file (optional, press Enter to skip): ").strip()
    label2 = None
    if fasta2:
        label2 = input("Enter label for second sequence (e.g., Horse): ").strip()

    # === Sequence 1 Analysis ===
    seq1 = read_fasta(fasta1)
    counts1 = count_codons(seq1)
    freqs1 = calculate_frequencies(counts1)
    export_to_csv(counts1, freqs1, label1, f"outputs/codon_usage_{label1}.csv")
    plot_codon_barplot(counts1, freqs1, label1, f"outputs/codon_plot_{label1}.png")
    plot_codon_heatmap(freqs1, label1, f"outputs/codon_heatmap_{label1}.png")

    # === Sequence 2 Analysis (if provided) ===
    if fasta2:
        seq2 = read_fasta(fasta2)
        counts2 = count_codons(seq2)
        freqs2 = calculate_frequencies(counts2)
        export_to_csv(counts2, freqs2, label2, f"outputs/codon_usage_{label2}.csv")
        plot_codon_barplot(counts2, freqs2, label2, f"outputs/codon_plot_{label2}.png")
        plot_codon_heatmap(freqs2, label2, f"outputs/codon_heatmap_{label2}.png")

        # === Comparison Outputs ===
        compare_codon_usage(freqs1, freqs2, label1, label2, "outputs/codon_comparison.png")
        codon_difference_heatmap(freqs1, freqs2, label1, label2, "outputs/codon_diff_heatmap.png")

    print("\n✅ Analysis complete! Results saved in 'outputs/' folder.")
