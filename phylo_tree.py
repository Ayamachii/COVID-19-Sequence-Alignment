import subprocess
from Bio import SeqIO, AlignIO
from tempfile import NamedTemporaryFile
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import csv


def align_sequences(case_sequences):
    """
    Align sequences using MAFFT via subprocess.

    Parameters:
    - case_sequences (list): List of SeqRecord objects representing sequences.

    Returns:
    - list: A list of aligned sequences as SeqRecord objects.
    """
    print("\nAligning sequences using MAFFT...")
    
    with NamedTemporaryFile(mode='w+', suffix=".fasta", delete=False) as tmp_input, \
         NamedTemporaryFile(mode='r+', suffix=".afa", delete=False) as tmp_output:
        
        # Write input sequences to a temporary FASTA file
        SeqIO.write(case_sequences, tmp_input.name, "fasta")
        tmp_input.flush()

        # Specify the full path to the MAFFT executable
        mafft_path = r"mafft-7.526-win64-signed\\mafft-win\\mafft.bat"
        
        # Construct the MAFFT command
        cmd = [mafft_path, "--auto", tmp_input.name]
        
        try:
            # Run the MAFFT command
            subprocess.run(cmd, check=True, stdout=tmp_output)
            tmp_output.flush()
            aligned_sequences = list(SeqIO.parse(tmp_output.name, "fasta"))
            print("Sequences aligned successfully.")
            return aligned_sequences
        
        except FileNotFoundError:
            print("MAFFT executable not found. Please ensure it is installed and in your PATH.")
            return []
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running MAFFT: {e}")
            return []

def print_distance_matrix(matrix):
    """
    Display a genetic distance matrix in the terminal in a formatted table.

    Parameters:
    - matrix: A Bio.Phylo.TreeConstruction.DistanceMatrix object.
    """
    print("\nDisplaying Genetic Distance Matrix in Terminal:")
    print("-" * 60)

    # Extract sequence names
    names = matrix.names
    name_width = max(len(name) for name in names) + 2  # Find the longest name for column width

    # Print header row
    header = f"{' ' * name_width}" + "".join(f"{name:>{name_width}}" for name in names)
    print(header)
    print("-" * len(header))

    # Print matrix rows with diagonal highlighted
    for i, (name, row) in enumerate(zip(names, matrix)):
        row_values = []
        for j, value in enumerate(row):
            if i == j:
                # Highlight diagonal values (self-comparisons)
                row_values.append(f"{value:.4f}")
                # row_values.append(f"\033[93m{value:.4f}\033[0m")  # Highlight in ipynb
            else:
                row_values.append(f"{value:.4f}")
        
        row_text = f"{name:<{name_width}}" + "".join(f"{val:>{name_width}}" for val in row_values)
        print(row_text)
    
    print("-" * 60)
    print("Genetic Distance Matrix displayed successfully.\n")


def list_clades(tree_file):
    tree = Phylo.read(tree_file, "newick")
    all_clades = list(tree.find_clades())
    print(f"Total number of clades: {len(all_clades)}")
    print("Clades details:")
    for i, clade in enumerate(all_clades, 1):
        descendants = [desc.name for desc in clade.get_terminals()]
        print(f"Clade {i}: Contains {len(descendants)} taxa")
        print("Descendants:", descendants)
    print("-" * 60)

def calculate_min_max_distances(matrix):
    dist_array = np.array(matrix)
    
    # Mask the diagonal for min calculations but reset for max calculations
    np.fill_diagonal(dist_array, np.inf)  # This is for ignoring during min calculations
    
    # Find the minimum distance
    min_distance = np.min(dist_array)
    
    # Get indices for the minimum distance
    min_indices = np.where(dist_array == min_distance)
    
    # Reset diagonal to a large number (not inf) for max calculation
    np.fill_diagonal(dist_array, -np.inf)
    
    # Find the maximum distance, ensuring to ignore the diagonal
    max_distance = np.max(dist_array[np.isfinite(dist_array)])  # Max, ignoring 'inf'
    
    # Get indices for the maximum distance
    max_indices = np.where(dist_array == max_distance)
    
    # Print results
    print("Minimum genetic distance between sequences:", min_distance)
    print(f"Closest sequences are between {matrix.names[min_indices[0][0]]} and {matrix.names[min_indices[1][0]]}.")
    print("Maximum genetic distance between sequences:", max_distance)
    print(f"Most distant sequences are between {matrix.names[max_indices[0][0]]} and {matrix.names[max_indices[1][0]]}.")
    print("-" * 60)




def create_combined_sequences(delta_selected, omicron_selected):
    """
    Joins the 10 Delta selected sequences with the selected Omicron sequences.

    Parameters:
    - delta_selected (list): List of selected Delta sequences.
    - omicron_selected (list): List of selected Omicron sequences.

    Returns:
    - list: Combined list of sequences.
    """
    print("\nCombining Delta and Omicron sequences...")
    combined = delta_selected + omicron_selected
    print("Sequences combined successfully.")
    return combined


def display_image(image_path):
    """
    Displays an image from a specified file path using matplotlib.

    Parameters:
    - image_path (str): Path to the image file.
    """
    print("\nDisplaying image...")
    img = mpimg.imread(image_path)
    plt.figure(figsize=(10, 8))
    plt.imshow(img)
    plt.axis('off')  # Turn off axis numbers and ticks
    plt.show()
    print("Image displayed successfully.\n")

def save_dict_to_csv(data, file_path):
    """
    Save a simple dictionary to a CSV file.

    Parameters:
        data (dict): Dictionary with key-value pairs.
        file_path (str): Path to save the CSV file.
    """
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Key', 'Value'])  # Header
        for key, value in data.items():
            writer.writerow([key, value])


if __name__ == "__main__":
    print("\nStarting Sequence Processing Pipeline\n")

    # Indices of the 10 chosen sequences from both files
    indices = [0, 1, 3, 6, 7, 8, 9, 11, 13, 14]
    
    print("Loading Delta and Omicron sequences...")
    delta_selected = [list(SeqIO.parse('Project_Files/Delta/gisaid_hcov-19_2021_12_31_12.fasta', "fasta"))[i] for i in indices]
    delta_cons = list(SeqIO.parse('Output_Files/delta_consensus.fasta', "fasta"))
    aligned_omicron = list(SeqIO.parse('Output_Files/aligned_omicron.fasta', "fasta"))
    omicron_sequences = list(SeqIO.parse('Project_Files/Omicron/gisaid_hcov-19_2021_12_31_12.fasta', "fasta"))
    omicron_selected = [omicron_sequences[i] for i in indices]
    print("Sequences loaded successfully.")

    # Rename the sequences for ease to see on the phylogeetic tree figure
    # use names_dict to map back to the original
    names_dict = {}
    for i, seq in enumerate(delta_selected):
        names_dict[f"Delta_{i+1}"] = seq.id
        seq.id = f"Delta_{i+1}"
        seq.name = seq.id
        seq.description = ""

    for i, seq in enumerate(omicron_selected):
        names_dict[f"Omicron_{i+1}"] = seq.id
        seq.id = f"Omicron_{i+1}"
        seq.name = seq.id
        seq.description = ""

    save_dict_to_csv(names_dict, 'Output_Files\\seq_names.csv')
    print("Sequence Names Dictionary saved as seq_names.csv")


    # Combine and align sequences
    combined_sequences = create_combined_sequences(delta_selected, omicron_selected)
    aligned_delta_omicron_sequences = align_sequences(combined_sequences)

    # Save aligned sequences
    output_file = "Output_Files\\aligned_all_delta_omicron_sequences.fasta"
    SeqIO.write(aligned_delta_omicron_sequences, output_file, "fasta")
    print(f"Aligned sequences saved to '{output_file}'.")

    # Calculate genetic distances
    print("\nCalculating genetic distances...")
    aligned_sequences_all = AlignIO.read(output_file, "fasta")
    distance_calculator = DistanceCalculator('identity') 
    distance_matrix = distance_calculator.get_distance(aligned_sequences_all)
    print("Genetic distance matrix calculated.")

    print_distance_matrix(distance_matrix)


    # Calculate minimum and maximum distances
    calculate_min_max_distances(distance_matrix)

    # Construct a phylogenetic tree
    print("\nConstructing phylogenetic tree using Neighbor-Joining method...")
    tree_constructor = DistanceTreeConstructor(distance_calculator, 'nj')
    phylogenetic_tree = tree_constructor.build_tree(aligned_sequences_all)
    print("Phylogenetic Tree Constructed.")

    print("\nListing clades from the constructed phylogenetic tree...")
    list_clades("Output_Files\\phylogenetic_tree.nwk")

    # Save and visualize the tree
    Phylo.write(phylogenetic_tree, "Output_Files\\phylogenetic_tree.nwk", "newick")
    print("Phylogenetic tree saved as 'phylogenetic_tree.nwk'.")

    print("\nDisplaying the Phylogenetic Tree...")
    fig, ax = plt.subplots(figsize=(14, 12))
    Phylo.draw(phylogenetic_tree, axes=ax, do_show=True)

    print("\nDisplaying the Phylogenetic Tree Using FigTree Graphical Viewer...")
    image_path = "Output_Files/Phylo_Tree.png"
    display_image(image_path)
    print("\nPipeline completed successfully")
