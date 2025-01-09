# 🧬 **Comparative Genomic Analysis of SARS-CoV-2 Delta and Omicron Variants in Ghana**

---

## 📚 **Project Overview**

This project focuses on a **comparative genomic analysis of SARS-CoV-2 Delta and Omicron variants in Ghana**, examining genetic differences, identifying dissimilar regions, and visualizing key findings using bioinformatics tools.

---

## 🛠️ **Setup Instructions**

### 1. **System Requirements**
- Python (≥3.8)  
- MAFFT version: **mafft-7.526-win64-signed**

### 2. **Environment Setup**
- Add the path of `mafft-7.526-win64-signed\mafft-win` to **System Variables**.

### 3. **Install Dependencies**
Run the following command in your project directory:
```bash
pip install -r requirements.txt
```

## 📂 **Output Files Explanation**

### 🧬 **FASTA Files**
- **delta_consesnsus.fasta:** A single record representing the **Delta consensus sequence** extracted from 10 sequences.  
- **aligned_omicron.fasta:** Contains **10 aligned Omicron sequences** from the MSA step.  
- **aligned_all_delta_omicron_sequences.fasta:** Alignment of **all 20 sequences** used for the **phylogenetic tree**.  
- **delta_omicron_aligned.fasta:** Alignment of **11 sequences** (Delta consensus + 10 Omicron sequences).  
- **dissimilar_cols.fasta:** Alignment of **dissimilar columns** extracted from `delta_omicron_aligned.fasta`.  
   - **Indices of these columns:** Found in `indices_dissimilar_cols.txt`.

---

### 📊 **CSV Files**
- **delta_group_analysis.csv:** Analysis results specific to **Delta regions**.  
- **omicron_group_analysis.csv:** Analysis results specific to **Omicron regions**.  
- **seq_names.csv:** Mapping of **sequence IDs** to their **names** in the phylogenetic tree.

---

### 📈 **PNG Files**
- **comparison_plot.png:** Visual comparison plot.  
- **insertion_nucleotide_frequencies.png:** Frequency of nucleotide **insertions**.  
- **deletion_nucleotide_frequencies.png:** Frequency of nucleotide **deletions**.  
- **grouped_vs_single_pie_chart.png:** Pie chart showing percentages of **dissimilar regions** vs **dissimilar residues**.  
- **mismatch_types_pie_chart.png:** Pie chart showing **deletions, insertions, and substitutions** percentages.  
- **Phylo_Tree.png:** Phylogenetic tree visualization.  
- **region_sizes_bar_plot.png:** Bar plot showing the **lengths of dissimilar regions**.  
- **substitution_heatmap.png:** Heatmap showing **substitution combinations percentages**.
- **single_columns_scatter_plot:** The scatter plot of the indices of SNPs in the alignment of 11 sequences.
- **regions_scatter_plot:** The scatter plot of the indices and lengths of dissimilar regions in the alignment of 11 sequences

---

### 📄 **Text Files**
- **indices_dissimilar_cols.txt:** Indices of **dissimilar columns** from `delta_omicron_aligned.fasta`.

---

### 📑 **Other Files**
- **phylo_tree.nwk.pdf:** Phylogenetic tree in **Newick format**.

---

## 🤝 **Contributors**
- **Aya Eyad**  
- **Shehab Mohamed**  
- **Farah Osama**  
- **Camellia Marwan**

**Supervisors:**  
- **Dr. Ibrahim Yousef**  
- **Eng. Yara Mahrous**

---
