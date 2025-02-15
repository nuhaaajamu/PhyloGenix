# PhyloGenix: Phylogenetic Tree Generator

PhyloGenix is a python based program that generates phylogenetic trees from sequence data in FASTA format. It leverages sequence alignment and evolutionary distance calculations to construct trees using UPGMA and Neighbor-Joining methods.
<br><br>

## Features
- Parses sequence data from FASTA files
- Computes pairwise evolutionary distances using sequence alignment
- Supports tree construction via UPGMA and Neighbor-Joining algorithms
- Provides ASCII and graphical visualizations of the phylogenetic tree
- Outputs the tree in Newick format for further analysis
<br><br>

## What is a FASTA File?
A FASTA file is a text-based format for representing nucleotide or protein sequences. Each sequence starts with a header line (beginning with >), followed by one or more lines of sequence data.
```plaintext
>sequence_1
ATCGATCGATCG
>sequence_2
GCTAGCTAGCTA
```
FASTA files are commonly used in bioinformatics for sequence alignment and phylogenetic analysis.
<br><br>

## Installation 
Ensure you have Python installed along with the required dependencies:
```plaintext
pip install biopython matplotlib
```

## How to Use
### 1. Running the Program
Run the script and provide user input when prompted:
```plaintext
python phlogenix.py
```
### 2. User Input
- Enter the path to your FASTA file
- Choose a tree construction method: UPGMA or NJ (Neighbor-Joining)
  
### Example Usage:
```plaintext
Enter the path to your FASTA file: example.fasta
Enter the tree construction method (UPGMA/NJ): NJ
```
### Output
- Displays the phylogenetic tree in ASCII format
- Outputs the tree in Newick format
- Visualizes the tree graphically
<br><br>

## Dependencies
- Python 3.7+
- Biopythonif 
- Matplotlib
<br><br>

## Future Improvements
- Support for additional phylogenetic algorithms
- Integration with multiple sequence alignment tools
- Enhanced interactive visualization options
<br><br>

## Contributions
Contributions are more than welcome! If this project is something that has sparked your interest, feel free to submit issues, feature requests, or pull requests.
<br><br>

## License
This is an open-source project and is licensed under the MIT License.






