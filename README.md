# Phylogenetic Tree Generator 🧬

PhyloGenix is a Python based program that generates phylogenetic trees from sequence data in FASTA format. It leverages sequence alignment and evolutionary distance calculations to construct trees using UPGMA and Neighbor-Joining algorithm methods.
<br><br>

## Features
- Parses sequence data from **FASTA** files
- Computes **pairwise evolutionary distances** using sequence alignment
- Supports tree construction via **UPGMA** and **Neighbor-Joining** algorithms
- Provides **ASCII** and **graphical** visualizations of the phylogenetic tree
- Outputs the tree in **Newick** format for further analysis
<p>&nbsp;</p> 

## Understanding the Methods Behind the Program

### What is a FASTA File?
A FASTA file is a text-based format for representing nucleotide or protein sequences. Each sequence starts with a header line 
**(beginning with >)**, followed by one or more lines of sequence data.

**Example Format:**
```plaintext
>sequence_1
ATCGATCGATCG
>sequence_2
GCTAGCTAGCTA
```
FASTA files are commonly used in bioinformatics for sequence alignment and phylogenetic analysis.
<p>&nbsp;</p> 

### UPGMA vs. Neighbor-Joining Algorithms

#### Key Differences:
| Feature            | Unwighted Pair Group Method with Arithmetic Mean (UPGMA)                              | Neighbor-Joining (NJ) |
|--------------------|----------------------------------|----------------------|
| **Methodology**   | Assumes a **molecular clock**, meaning all lineages evolve at a constant rate. | Does **not** assume a molecular clock, allowing for varying rates of evolution. |
| **Tree Type**     | Produces an **ultrametric tree**, where all leaf nodes are equidistant from the root. | Produces an **additive tree**, where branch lengths reflect evolutionary distance. |
| **Accuracy**      | Works best for **data with equal mutation rates** across lineages. | More flexible and accurate for **real-world biological data** with different mutation rates. |
| **Visualization** | The tree appears more **balanced and uniform** due to the equal evolutionary rates assumption. | The tree may appear **asymmetric**, reflecting varying rates of evolution across species. |
<p>&nbsp;</p>

#### How This Affects Visualization:
- **UPGMA Trees**: Tend to look more **symmetrical**, as all taxa are assumed to have evolved at the same rate.
- **NJ Trees**: Can appear **asymmetrical**, as taxa with different mutation rates will have branches of varying lengths.
<p>&nbsp;</p> 

### Newick (NHX) Format  

The **Newick format** is a standard text-based notation for representing phylogenetic trees. It expresses evolutionary relationships using a **nested, parenthetical format** where taxa are enclosed in parentheses and branch lengths are represented as numbers.  

#### **Example Format:**  
```plaintext
(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
```
In this example:
- A and B share a common ancestor, with branch lengths 0.1 and 0.2.
- C and D share a common ancestor, with branch lengths 0.3 and 0.4.
- The root node connects all branches at a distance of 0.5.
<p>&nbsp;</p> 
Newick format is widely used in phylogenetics for storing and exchanging tree structures across different bioinformatics tools.
<p>&nbsp;</p> 

## How to Run the Program
### 1. Clone the Repository
```plaintext
git clone https://github.com/nuhaaajamu/phylogenix.git
cd phylogenix
```
### 2. Install Dependencies
Ensure you have Python (version 3.7+) installed, then install these required libraries:
```plaintext
pip install biopython matplotlib
```
### 3. Run Program
Run the script and provide user input when prompted:
```plaintext
python phylo_tree_builder.py
```
### 4. User Input
- Enter the path to your FASTA file
- Choose a tree construction method: UPGMA or NJ (Neighbor-Joining)

### Example Usage:
```plaintext
Enter the path to your FASTA file: file_name.fasta
Select a tree construction algorithm (UPGMA/NJ): UPGMA
```
<p>&nbsp;</p> 

## Visualization Using FASTA File Example

Example trees generated using the **example_file.fasta** file in the repository 
<br><br>

<p align="center"><strong>UPGMA Algorithm Tree Construction</strong></p>

<p align="center">
  <img src="example_tree_UPGMA.png" alt="Generated Phylogenetic Tree" width="60%">
</p>

<br><br>

<p align="center"><strong>Neighbor-Joining Algorithm Tree Construction</strong></p>

<p align="center">
  <img src="example_tree_NJ.png" alt="Generated Phylogenetic Tree" width="60%">
</p>
<p>&nbsp;</p> 

## Dependencies
- Python 3.7+
- Biopython 
- Matplotlib

## Future Improvements
- Support for additional phylogenetic algorithms
- Integration with multiple sequence alignment tools
- Enhanced interactive visualization options

## Contributions
Contributions are more than welcome! If this project is something that has sparked your interest, feel free to submit issues, feature requests, or pull requests.

## License
This is an open-source project and is licensed under the MIT License.






