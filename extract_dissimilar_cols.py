import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import matplotlib.pyplot as plt

def analyze_and_visualize_mismatches(alignment, mode_mismatch_positions):
    """
    Analyze mismatched columns and classify them into Insertion, Deletion, and Substitution.
    Create a pie chart showing the percentage of each mismatch type.

    Parameters:
    - alignment (MultipleSeqAlignment): The alignment object.
    - mode_mismatch_positions (list): Indices of mismatched columns.

    Saves:
    - A pie chart showing mismatch types as 'mismatch_types_pie_chart.png'.
    """
    print("\nAnalyzing mismatched columns for Insertion, Deletion, and Substitution...")

    delta_consensus = alignment[-1].seq  # Last sequence is Delta consensus
    omicron_sequences = alignment[:-1]  # Other sequences are Omicron

    # Initialize mismatch counters
    mismatch_counts = {"Insertion": 0, "Deletion": 0, "Substitution": 0}

    for i in mode_mismatch_positions:
        delta_char = delta_consensus[i]
        omicron_column = [record.seq[i] for record in omicron_sequences]
        mode_char, _ = Counter(omicron_column).most_common(1)[0]

        if mode_char == '-' and delta_char != '-':
            mismatch_counts["Insertion"] += 1
        elif mode_char != '-' and delta_char == '-':
            mismatch_counts["Deletion"] += 1
        elif mode_char != '-' and delta_char != '-' and mode_char != delta_char:
            mismatch_counts["Substitution"] += 1

    print("\nMismatch Analysis Complete:")
    for mismatch_type, count in mismatch_counts.items():
        print(f"   - {mismatch_type}: {count} columns")

    # Calculate percentages
    total_mismatches = sum(mismatch_counts.values())
    mismatch_percentages = {k: (v / total_mismatches) * 100 for k, v in mismatch_counts.items()}

    # Plot the pie chart
    print("\nGenerating mismatch type pie chart...")
    labels = list(mismatch_percentages.keys())
    sizes = list(mismatch_percentages.values())
    colors = ['skyblue', 'lightcoral', 'gold']

    plt.figure(figsize=(8, 8))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
    plt.title('Distribution of Mismatch Types (Insertion, Deletion, Substitution)')
    plt.savefig('Output_Files\\mismatch_types_pie_chart.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Pie chart saved to 'Output_Files\\mismatch_types_pie_chart.png'.")


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
    output_fasta_file='Output_Files\\dissimilar_cols.fasta'
):
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

    print("\n🔹 Analyzing mismatch types...")
    analyze_and_visualize_mismatches(alignment, consensus_indices)

    print('\n')
    _, _ = extract_columns_consensus_vs_mode(alignment, exclude_unknowns=False)

    print("\nSequence processing pipeline completed successfully!\n")
