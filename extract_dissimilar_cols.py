import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import Counter, defaultdict
import pandas as pd
import random

from collections import defaultdict, Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def analyze_and_visualize_mismatches(alignment, mode_mismatch_positions, save_path='Output_Files'):
    """
    Analyze mismatched columns and classify them into Insertion, Deletion, and Substitution.
    Calculate nucleotide frequencies for Insertions and Deletions.
    Generate a heatmap for Substitutions and visualize all results.
    For Justifying getting the mode (as if getting the Omicron consensus) to compare to the delta consesnsus
    This is more logical than comparing if there was any discrepancy between any 2 of the 11 characters
    since if there was a descrepancy between two of the omicron sequences, it dorsn't indicate a difference between
    omicron and delta, which is our interest here.

    Parameters:
    - alignment (MultipleSeqAlignment): The alignment object.
    - mode_mismatch_positions (list): Indices of mismatched columns.
    - save_path (str): Directory path to save the plots.

    Saves:
    - Pie chart for mismatch types.
    - Bar plots for insertion and deletion nucleotide frequencies.
    - Heatmap for substitution frequencies.
    """
    print("\nAnalyzing mismatched columns for Insertion, Deletion, and Substitution...")

    delta_consensus = alignment[-1].seq  # Last sequence is Delta consensus
    omicron_sequences = alignment[:-1]  # Other sequences are Omicron

    # Initialize mismatch counters
    mismatch_counts = {"Insertion": 0, "Deletion": 0, "Substitution": 0}
    insertion_nucleotides = Counter()
    deletion_nucleotides = Counter()
    substitution_counts = defaultdict(int)

    for i in mode_mismatch_positions:
        # get chracter at position i in delta consensus
        delta_char = delta_consensus[i]

        # get the most common character in the 10 aligned Omicron seqs (as if consensus)
        omicron_column = [record.seq[i] for record in omicron_sequences]
        mode_char, _ = Counter(omicron_column).most_common(1)[0]

        # Count each case for each nucleotide
        if mode_char == '-' and delta_char != '-':
            mismatch_counts["Insertion"] += 1
            insertion_nucleotides[delta_char] += 1
        elif mode_char != '-' and delta_char == '-':
            mismatch_counts["Deletion"] += 1
            deletion_nucleotides[mode_char] += 1
        elif mode_char not in ['-', 'n', 'N'] and delta_char not in ['-', 'n', 'N'] and mode_char != delta_char:
            mismatch_counts["Substitution"] += 1
            substitution_counts[(delta_char, mode_char)] += 1

    print("\nMismatch Analysis Complete:")
    for mismatch_type, count in mismatch_counts.items():
        print(f"   - {mismatch_type}: {count} columns")

    # Calculate percentages
    total_mismatches = sum(mismatch_counts.values())
    mismatch_percentages = {k: (v / total_mismatches) * 100 for k, v in mismatch_counts.items()}

    # Pie Chart for Mismatch Types
    print("\nGenerating mismatch type pie chart...")
    labels = list(mismatch_percentages.keys())
    sizes = list(mismatch_percentages.values())
    colors = ['skyblue', 'lightcoral', 'gold']

    plt.figure(figsize=(8, 8))
    _, _, autotexts = plt.pie(
        sizes, 
        labels=labels, 
        autopct='%1.1f%%', 
        startangle=140, 
        colors=colors,
        textprops={'fontsize': 14}
    )
    for autotext in autotexts:
        autotext.set_fontsize(18)
        autotext.set_color('black')

    plt.title('Distribution of Mismatch Types (Insertion, Deletion, Substitution)')
    plt.savefig(os.path.join(save_path, 'mismatch_types_pie_chart.png'), dpi=300, bbox_inches='tight')
    plt.show()

    # Bar Charts for Insertions and Deletions of different nucleotides
    print("\nGenerating insertion and deletion bar charts...")
    for label, nucleotide_counts in zip(["Insertion", "Deletion"], [insertion_nucleotides, deletion_nucleotides]):
        plt.figure(figsize=(8, 6))
        nucleotides, counts = zip(*nucleotide_counts.items()) if nucleotide_counts else ([], [])
        plt.bar(nucleotides, counts)
        plt.title(f'{label} Nucleotide Frequencies')
        plt.xlabel('Nucleotide', fontsize=14)
        plt.ylabel('Count', fontsize=14)
        plt.savefig(os.path.join(save_path, f'{label.lower()}_nucleotide_frequencies.png'), dpi=300, bbox_inches='tight')
        plt.show()

    # Print raw substitution counts
    print("\nSubstitution Counts Breakdown:")
    total_substitutions = sum(substitution_counts.values())
    if total_substitutions == 0:
        print("\nNo substitutions detected! Heatmap will remain empty.")
        return

    # print the counts for every combination in the terminal
    for (sub, orig), count in substitution_counts.items():
        print(f"{orig} → {sub}: {count}")

    # Create substitution matrix with percentages
    substitution_matrix = pd.DataFrame(0, index=nucleotides, columns=nucleotides, dtype=float)
    for (sub, orig), count in substitution_counts.items():
        if orig in nucleotides and sub in nucleotides:
            substitution_matrix.at[orig, sub] = (count / total_substitutions) * 100

    # Print substitution matrix
    print("\nSubstitution Matrix (Percentages):")
    print(substitution_matrix)

    # Plot heatmap
    sns.heatmap(
        substitution_matrix,
        annot=True,
        cmap='Blues',
        fmt='.2f',
        cbar_kws={'label': 'Percentage (%)'},
        annot_kws={"size": 14}  # Increase font size for annotations
    )
    plt.title('Substitution Frequency Heatmap (Percentages)', fontsize=16)
    plt.xlabel('Substituted Nucleotide', fontsize=14)
    plt.ylabel('Original Nucleotide', fontsize=14)

    # Save Heatmap to png file
    os.makedirs(save_path, exist_ok=True)
    plt.savefig(os.path.join(save_path, 'substitution_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.show()

    print(f"\nHeatmap saved as '{os.path.join(save_path, 'substitution_heatmap.png')}'")

    print(f"\nAll visualizations saved in '{save_path}'. Analysis complete!")


def extract_columns_consensus_vs_mode(alignment, exclude_unknowns=False):
    """
    Extract columns where the consensus sequence differs from the mode character
    in the Omicron sequences and return column indices. (As if comparing the 2 consensuses of Delta and Omicron)
    For the exclude_unknowns param: We think excluding where the consensus has N is more logical since it does not mean a
    descrepancy between the 2 variants but rather missing information in the delta sewuences. 
    Also notice how they are about 600 columns vs actual variants which are about 300. (Will be printed afterwards)

    Parameters:
    - alignment (MultipleSeqAlignment): Aligned sequences, with the consensus sequence as the last sequence.
    - exclude_unknowns (bool): Whether to exclude columns where the Delta consensus contains 'N' or 'n'. Deafult is false

    Returns:
    - tuple: A tuple containing:
        - List of column indices with mismatches.
        - List of tuples with (record_id, extracted_column_sequence).
    """
    print(f"Extracting mismatched columns (exclude_unknowns={exclude_unknowns})...")

    consensus_sequence = alignment[-1].seq  # Last sequence is the delta consensus
    omicron_sequences = alignment[:-1]  # Remaining sequences are Omicron

    mode_mismatch_positions = []  # Store column indices with consensus-mode mismatches to be exported to a txt file

    for i in range(alignment.get_alignment_length()):
        if exclude_unknowns and consensus_sequence[i] in ('N', 'n'):
            continue

        # Get the most common (mode) char in the aligned omicron to compare to the consensus char at the same location
        omicron_column = [record.seq[i] for record in omicron_sequences]
        mode_char, _ = Counter(omicron_column).most_common(1)[0]

        if consensus_sequence[i] != mode_char:
            mode_mismatch_positions.append(i)

    # Extract these dissimilar columns for each sequence
    extracted_columns = [
        (record.id, ''.join(record.seq[i] for i in mode_mismatch_positions))
        for record in alignment
    ]

    print(f"Found {len(mode_mismatch_positions)} mismatched columns.")
    return mode_mismatch_positions, extracted_columns

def save_to_fasta_with_indices(
    extracted_columns, column_indices, 
    index_file='Output_Files\\indices_dissimilar_cols.txt', 
    output_fasta_file='Output_Files\\dissimilar_cols.fasta'):
    """
    Save extracted columns to a fasta file and their indices in the 11 seq alignment to a separate txt file.

    Parameters:
    - extracted_columns (list): List of tuples (record_id, sequence).
    - column_indices (list): List of column indices extracted.
    - index_file (str): File path to save indices (.txt).
    - output_fasta_file (str): File path to save extracted columns (.fasta).
    """
    print(f"Saving column indices to '{index_file}'...")
    with open(index_file, 'w') as index_f:
        index_f.write("# Extracted Column Indices\n")
        # Write the indeces successively to the text file on one line
        index_f.write(', '.join(map(str, column_indices)) + '\n')
    print(f"Column indices saved to '{index_file}'.")

    print(f"Saving extracted columns to '{output_fasta_file}'...")
    records = [
        SeqRecord(Seq(seq), id=record_id, description="Extracted Columns")
        for record_id, seq in extracted_columns
    ]
    SeqIO.write(records, output_fasta_file, "fasta")
    print(f"Extracted columns saved to '{output_fasta_file}'.")

def run_mafft_add(delta_file, omicron_file, output_file):
    """
    Run MAFFT to add a single consensus sequence to an existing alignment.

    Parameters:
    - delta_file (str): Path to the consensus Delta sequence file.
    - omicron_file (str): Path to the aligned Omicron sequences file.
    - output_file (str): Path to save the final alignment.

    Returns:
    - list: A list of aligned sequences as SeqRecord objects.
    """
    print(f"\nRunning MAFFT to add '{delta_file}' to '{omicron_file}'...")

    mafft_path = r"mafft-7.526-win64-signed\\mafft-win\\mafft.bat"

    # This add flag adds a single sequence to an existing alignment without modifying the existing alignment
    cmd = [mafft_path, "--add", delta_file, omicron_file]

    try:
        with open(output_file, "w") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)
        print(f"MAFFT alignment completed. Output saved to '{output_file}'.")

        aligned_sequences = list(SeqIO.parse(output_file, "fasta"))
        return aligned_sequences

    except FileNotFoundError:
        print("MAFFT executable not found. Please ensure it is installed and in your PATH.")
        return []
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running MAFFT: {e.stderr}")
        return []

def group_indices(indices):
    """
    Group the indices of dissimilar columns into contiguous regions and separate single columns.

    Parameters:
    - indices (list): List of indices of dissimilar columns separately.

    Returns:
    - regions (list): Lists of indices, each representing one contiguous region.
    - len(regions) (int): The number of regions found (more than 2 successive columns).
    - single_columns_count (int): The number of single columns found.
    - single_columns_indices (list): Indices of single columns.
    """
    if not indices:
        return [], 0, 0, []
    
    indices = sorted(set(indices))  # Ensure sorted and no duplicates
    regions = []  # List to store contiguous regions
    single_columns_count = 0  # Count of single columns
    single_columns_indices = []  # List to store indices of single columns
    
    current_region = [indices[0]]
    
    for i in range(1, len(indices)):
        if indices[i] == indices[i-1] + 1:
            # Successive index, add to current region
            current_region.append(indices[i])
        else:
            # End of a region
            if len(current_region) == 1:
                single_columns_count += 1
                single_columns_indices.append(current_region[0])
            else:
                regions.append(current_region)
            current_region = [indices[i]]
    
    # Handle the last region similarly
    if len(current_region) == 1:
        single_columns_count += 1
        single_columns_indices.append(current_region[0])
    else:
        regions.append(current_region)
    
    return regions, len(regions), single_columns_count, single_columns_indices

def visualize_regions(regions, num_regions, single_columns_count, single_columns_indeces, save_path='Output_Files\\'):
    """
    Visualizes statistics and insights about grouped regions and single columns.
    
    Parameters:
        regions (list): A list of grouped regions (each region is a list of indices).
        num_regions (int): The number of grouped regions.
        single_columns (int): The number of single isolated indices of dissimilar sinle columns (SNPs).
        save_path (str): The directory path where figures will be saved. Default is 'figures'.
    
    """
    # Statistics
    region_sizes = [len(region) for region in regions]

    max_region_size = max(region_sizes) if region_sizes else 0
    min_region_size = min(region_sizes) if region_sizes else 0
    total_count = sum(region_sizes) + single_columns_count
    coverage_by_regions = sum(region_sizes)
    percentage_single_columns = (single_columns_count / total_count) * 100 if total_count > 0 else 0
    
    # Print to termial
    print(f"Number of Originally Extracted Indices: {total_count}")
    print(f"Number of Regions: {num_regions}")
    print(f"Number of Single Columns: {single_columns_count}")
    print(f"Largest Region Size: {max_region_size}")
    print(f"Smallest Region Size: {min_region_size}")
    print(f"Coverage by Regions: {coverage_by_regions}")
    print(f"Percentage of Single Columns: {percentage_single_columns:.2f}%")
    
    # Bar Plot of Region Sizes
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(region_sizes)), region_sizes)
    plt.title('Sizes of Grouped Regions')
    plt.xlabel('Region Index')
    plt.ylabel('Size of Region')
    plt.savefig(os.path.join(save_path, 'region_sizes_bar_plot.png'))
    plt.show()
    
    # Pie Chart of the length covered by Single Columns vs Grouped Regions
    labels = ['Grouped Indices', 'Single Columns']
    sizes = [coverage_by_regions, single_columns_count]
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%')
    plt.title('Proportion of Grouped Indices vs Single Columns')
    plt.savefig(os.path.join(save_path, 'grouped_vs_single_pie_chart.png'))
    plt.show()
    
    # Scatter Plot Single Columns
    plt.figure(figsize=(10, 3))
    plt.scatter(single_columns_indeces, [1] * len(single_columns_indeces), 
                label='Single Columns', 
                marker='x', 
                color='red', 
                s=50)
    plt.title('Single Columns Distribution')
    plt.xlabel('Index Value')
    plt.yticks([])

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), frameon=False)

    # Save and Display Single Columns Plot
    plt.savefig(os.path.join(save_path, 'single_columns_scatter_plot.png'), dpi=300, bbox_inches='tight')
    plt.show()


    # Scatter Plot of Regions with Different Colors
    plt.figure(figsize=(10, 3))
    # Generate random distinct colors for each region
    region_colors = [f'#{random.randint(0, 0xFFFFFF):06x}' for _ in regions]

    # Plot each region and add to legend with start index and length
    legend_handles = []
    for i, region in enumerate(regions):
        plt.scatter(region, [1] * len(region), 
                    color=region_colors[i], 
                    marker='o', 
                    s=50)
        legend_label = f'Region {region[0]} (Length: {len(region)})'
        legend_handles.append(plt.Line2D([0], [0], marker='o', color=region_colors[i], lw=0, label=legend_label))

    plt.title('Regions Distribution')
    plt.xlabel('Index Value')
    plt.yticks([])

    # Custom Legend to carry the starting index and length of every region
    plt.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)

    # Save and Display Regions Plot
    plt.savefig(os.path.join(save_path, 'regions_scatter_plot.png'), dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    print("\nStarting Sequence Processing Pipeline\n")

    # Predefined indices for sequence selection
    indices = [0, 1, 3, 6, 7, 8, 9, 11, 13, 14]
    
    delta_selected = [list(SeqIO.parse('Project_Files/Delta/gisaid_hcov-19_2021_12_31_12.fasta', "fasta"))[i] for i in indices]
    delta_cons = list(SeqIO.parse('Output_Files/delta_consensus.fasta', "fasta"))
    aligned_omicron = list(SeqIO.parse('Output_Files/aligned_omicron.fasta', "fasta"))

    print(f"Total Delta sequences selected: {len(delta_selected)}")
    print(f"Total Omicron sequences aligned: {len(aligned_omicron)}")

    run_mafft_add(
        delta_file="Output_Files\\delta_consensus.fasta",
        omicron_file="Output_Files\\aligned_omicron.fasta",
        output_file="Output_Files\\delta_omicron_aligned.fasta"
    )

    alignment = AlignIO.read("Output_Files\\delta_omicron_aligned.fasta", "fasta")

    print('\n')
    consensus_indices, extracted = extract_columns_consensus_vs_mode(alignment, exclude_unknowns=True)
    print('\n')
    save_to_fasta_with_indices(extracted, consensus_indices)

    print("\nAnalyzing mismatch types...")
    analyze_and_visualize_mismatches(alignment, consensus_indices)

    # Group Indices to regions, single columns and count the number of each
    regions, num_regions, single_columns_count, single_columns_indeces = group_indices(consensus_indices)
    
    visualize_regions(regions, num_regions, single_columns_count, single_columns_indeces)

    print('\n')
    _, _ = extract_columns_consensus_vs_mode(alignment, exclude_unknowns=False)

    print("\nSequence processing pipeline completed successfully!\n")
