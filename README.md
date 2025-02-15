# PhyloGenix: Phylogenetic Tree Generator

PhyloGenix is a python based program that generates phylogenetic trees from sequence data in FASTA format. It leverages sequence alignment and evolutionary distance calculations to construct trees using UPGMA and Neighbor-Joining methods.
<br><br>

## Features
- Parses sequence data from FASTA files
- Computes pairwise evolutionary distances using sequence alignment
- Supports tree construction via UPGMA and Neighbor-Joining algorithms
- Provides ASCII and graphical visualizations of the phylogenetic tree
- Outputs the tree in Newick format for further analysis

## What is a FASTA File?
A FASTA file is a text-based format for representing nucleotide or protein sequences. Each sequence starts with a header line 
**(beginning with >)**, followed by one or more lines of sequence data.
```plaintext
>sequence_1
ATCGATCGATCG
>sequence_2
GCTAGCTAGCTA
```
FASTA files are commonly used in bioinformatics for sequence alignment and phylogenetic analysis.

## Installation 
Ensure you have Python installed along with the required dependencies:
```plaintext
pip install biopython matplotlib
```

## How to Use
### 1. Running the Program
Run the script and provide user input when prompted:
```plaintext
python phylo_tree_builder.py
```
### 2. User Input
- Enter the path to your FASTA file
- Choose a tree construction method: UPGMA or NJ (Neighbor-Joining)
  
### Example Usage:
```plaintext
Enter the path to your FASTA file: example.fasta
Enter the tree construction method (UPGMA/NJ): UPGMA
```
### Example Output

Tree visualization using the provided "example_fasta.fasta" file:

![Generated Phylogenetic Tree](example_tree.png)

To see this output yourself, run:
```bash
python phylo_tree_builder.py example_fasta.fasta
```

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






