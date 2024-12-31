from Bio import SeqIO
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
import os


def analyze_nucleotide_content(fasta_file):
    """
    Analyze nucleotide content and percentages for sequences in a FASTA file.

    parameters:
        fasta_file (str): Path to the FASTA file containing sequences.

    Returns:
        pd.DataFrame: DataFrame with nucleotide counts, percentages, and CG content for each sequence,
                    and overall averages for the group.
    """
    # List to hold individual sequence data
    data = []
    # Process each sequence in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        length = len(sequence)

        # Count nucleotides
        count_a = sequence.count("A")
        count_t = sequence.count("T")
        count_g = sequence.count("G")
        count_c = sequence.count("C")

        # Calculate nucleotide percentagesS
        percent_a = (count_a / length) * 100
        percent_t = (count_t / length) * 100
        percent_g = (count_g / length) * 100
        percent_c = (count_c / length) * 100

        # Calculate CG content
        cg_content = percent_c + percent_g

        # Append data for this sequence
        data.append({
            "Sequence_ID": record.id,
            "A_percentage": percent_a,
            "T_percentage": percent_t,
            "G_percentage": percent_g,
            "C_percentage": percent_c,
            "CG_content": cg_content
        })

    # Create DataFrame for individual sequence data
    df = pd.DataFrame(data)

    # Calculate averages across all sequences
    avg_data = {
        "Sequence_ID": "Average",
        "A_percentage": df["A_percentage"].mean(),
        "T_percentage": df["T_percentage"].mean(),
        "G_percentage": df["G_percentage"].mean(),
        "C_percentage": df["C_percentage"].mean(),
        "CG_content": df["CG_content"].mean()
    }
    # Append average row to DataFrame
    df = pd.concat([df, pd.DataFrame([avg_data])], ignore_index=True)
    return df


delta_file = "./Project_Files/Delta/gisaid_hcov-19_2021_12_31_12.fasta"
omicron_file = "./Project_Files/Omicron/gisaid_hcov-19_2021_12_31_12.fasta"

df_delta = analyze_nucleotide_content(delta_file)
df_omicron = analyze_nucleotide_content(omicron_file)

# Save results to a CSV file 
# df_delta.to_csv("./Output_Files/delta_group_analysis.csv", index=False)
# df_omicron.to_csv("./Output_Files/omicron_group_analysis.csv", index=False)


omicron_avg = df_omicron.iloc[-1].drop("Sequence_ID")
delta_avg = df_delta.iloc[-1].drop("Sequence_ID")
result_of_comparison_df=pd.DataFrame({
        "Group": ["Delta", "Omicron"],
        "A_percentage": [delta_avg["A_percentage"], omicron_avg["A_percentage"]],
        "T_percentage": [delta_avg["T_percentage"], omicron_avg["T_percentage"]],
        "G_percentage": [delta_avg["G_percentage"], omicron_avg["G_percentage"]],
        "C_percentage": [delta_avg["C_percentage"], omicron_avg["C_percentage"]],
        "CG_content": [delta_avg["CG_content"], omicron_avg["CG_content"]]
    })


# print(tabulate(df_results_omicron, headers='keys', tablefmt='pretty', showindex=False))
# print(tabulate(df_results_delta, headers='keys', tablefmt='pretty', showindex=False))
print(tabulate(result_of_comparison_df, headers='keys', tablefmt='pretty', showindex=False))


###########################################################################################################


# Extract data for visualization
groups = result_of_comparison_df["Group"]
categories = ["A_percentage", "T_percentage", "G_percentage", "C_percentage", "CG_content"]

# Data for bar chart
data = result_of_comparison_df[categories].values.T  # Transpose for grouped bar chart
x = np.arange(len(categories))  # Category positions

# Plotting grouped bar chart
fig, ax = plt.subplots(figsize=(10, 6))
width = 0.35  # Bar width

# Bars for each group
ax.bar(x - width/2, data[:, 0], width, label=groups[0])
ax.bar(x + width/2, data[:, 1], width, label=groups[1])

# Customization
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylabel("Percentage (%)")
ax.set_title("Average Nucleotide Percentages and CG Content by Group")
ax.legend(title="Group")
#save the image
output_path = "./Output_Files/comparison_plot.png"  
plt.savefig(output_path)
# Show plot
plt.tight_layout()
plt.show()
