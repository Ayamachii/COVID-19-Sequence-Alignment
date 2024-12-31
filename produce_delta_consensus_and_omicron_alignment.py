import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tempfile import NamedTemporaryFile
from collections import Counter


def get_consensus(ref_sequences):
    """
    Construct a consensus sequence from a list of reference sequences.

    Parameters:
    - ref_sequences (list): A list of reference sequences (strings).

    Returns:
    - str: A consensus sequence constructed by taking the most frequent character
           at each position, ignoring gaps ('-').
    """
    consensus = []
    for i in range(len(ref_sequences[0])):
        # Collect the characters at position i, ignoring gaps
        chars = [seq[i] for seq in ref_sequences if seq[i] != '-']
        if chars:
            # Choose the most common character (mode)
            consensus.append(max(set(chars), key=chars.count))
        else:
            # If only gaps are present, keep the gap in the consensus
            consensus.append('-')
    return ''.join(consensus)


def align_sequences(case_sequences):
    """
    Align sequences using MAFFT, used to align sequences not a sequence and an alignment.

    Parameters:
    - case_sequences (list): List of SeqRecord objects of sequences to align.

    Returns:
    - list: A list of aligned sequences as SeqRecord objects.
    """
    with NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as tmp_input, \
         NamedTemporaryFile(mode='r+', suffix=".afa", delete=False) as tmp_output:
        
        # Write input sequences to a temporary FASTA file
        SeqIO.write(case_sequences, tmp_input.name, "fasta")
        tmp_input.flush()

        # Specify the path to the MAFFT executable
        mafft_path = r"mafft-7.526-win64-signed\\mafft-win\\mafft.bat"
        
        # The mafft command to run
        cmd = [mafft_path, "--auto", tmp_input.name]
        
        try:
            # Run the MAFFT alignment using subprocess to run the command 
            subprocess.run(cmd, check=True, stdout=tmp_output)
            tmp_output.flush()
            aligned_sequences = list(SeqIO.parse(tmp_output.name, "fasta"))
            return aligned_sequences
        
        except FileNotFoundError:
            print("MAFFT executable not found. Please ensure it is installed and in your PATH.")
            return []
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running MAFFT: {e}")
            return []


def alignment_stats(alignment):
    """
    Calculate and print statistics for a given alignment.

    Parameters:
    - alignment (list): A list of SeqRecord objects representing aligned sequences.

    Prints:
    - Average match rate across alignment columns.
    """
    # Transpose alignment columns to loop on cols instead of sequneces
    alignment_columns = zip(*[str(seq.seq) for seq in alignment])
    
    # Calculate match rate for each column
    column_match_rate = [
        max(Counter(col).values()) / len(col) for col in alignment_columns
    ]
    
    # Calculate average match rate on all columns
    avg_match_rate = sum(column_match_rate) / len(column_match_rate)
    print(f"Average match rate across alignment columns: {avg_match_rate:.2f}")


def process_files(delta_path, omicron_path):
    """
    Process Delta and Omicron variant sequences:
    - Extract specific sequences.
    - Generate consensus sequence for Delta.
    - Align Omicron sequences.
    - Calculate alignment statistics.

    Parameters:
    - delta_path (str): Path to the Delta variant sequence file.
    - omicron_path (str): Path to the Omicron variant sequence file.
    """
    # Load sequences from the provided files
    delta_sequences = list(SeqIO.parse(delta_path, "fasta"))
    omicron_sequences = list(SeqIO.parse(omicron_path, "fasta"))

    # Select 10 sequences based on predefined indices
    delta_selected = [delta_sequences[i] for i in indices]
    omicron_selected = [omicron_sequences[i] for i in indices]

    print('*'*100)
    print('\n')
    print('*'*100)

    print("Aligning Omicron Sequences...")
    # Align Omicron sequences
    aligned_omicron = align_sequences(omicron_selected)
    print("Alignment of Omicron Done")

    # Save aligned Omicron sequences to a file
    if aligned_omicron:
        with open("Output_Files\\aligned_omicron.fasta", "w") as output_file:
            SeqIO.write(aligned_omicron, output_file, "fasta")
        print("Aligned Omicron sequences saved to aligned_omicron.fasta")
    else:
        print("No aligned sequences were generated.")
    
    print('*'*50)   
    print('\n')
    
    # Display alignment statistics
    print("Omicron Alignment Statistics:")
    alignment_stats(aligned_omicron)

    print('*'*100)
    print('\n')
    print('*'*100)

    print('Getting Delta Consensus...')
    # Construct Delta consensus sequence
    delta_consensus = get_consensus([str(seq.seq) for seq in delta_selected])

    # Create a SeqRecord object for the consensus sequence
    delta_consensus_record = SeqRecord(
        Seq(delta_consensus),
        id="delta_consensus",
        description="delta_consensus_Ghana"
    )

    # Save the Delta consensus sequence to a file
    with open("Output_Files\\delta_consensus.fasta", "w") as output_file:
        SeqIO.write(delta_consensus_record, output_file, "fasta")
    print("Delta Consensus Sequence saved to delta_consensus.fasta")
    print('*'*100)


# -------------------------------
# Main Function
# -------------------------------
if __name__ == "__main__":
    # Predefined indices for sequence selection
    indices = [0, 1, 3, 6, 7, 8, 9, 11, 13, 14]
    
    # Define paths to input files
    delta_path = 'Project_Files/Delta/gisaid_hcov-19_2021_12_31_12.fasta'
    omicron_path = 'Project_Files/Omicron/gisaid_hcov-19_2021_12_31_12.fasta'
    
    # Process Delta and Omicron variant sequences: Extract specific sequences, Generate consensus sequence for Delta.
    # Align Omicron sequences, Calculate alignment statistics.
    process_files(delta_path, omicron_path)
