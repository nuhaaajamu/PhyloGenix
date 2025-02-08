from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
import matplotlib.pyplot as plt


def parse_fasta_file(file_path):
    parsed_sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            # Retrieve the sequence identifier and corresponding sequence.
            parsed_sequences[record.id] = str(record.seq)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        exit(1)
    except Exception as e:
        print(f"Error reading file '{file_path}': {e}")
        exit(1)

    if not parsed_sequences:
        print(f"Error: The file '{file_path}' is empty or not in valid FASTA format.")
        exit(1)

    return parsed_sequences


def calculate_distance_matrix(parsed_sequences):
    seq_aligner = PairwiseAligner()
    seq_ids = list(parsed_sequences.keys())
    id_num = len(seq_ids)
    lower_triangular_matrix = []

    # Generate a lower triangular matrix of pairwise distances for compatibility with DistanceTreeConstructor.
    for i in range(id_num):
        row = []
        for j in range (i + 1):
            if i == j:
                # The distance of a sequence to itself is always zero.
                row.append(0.0)
            else:
                seq_one = parsed_sequences[seq_ids[i]]
                seq_two = parsed_sequences[seq_ids[j]]

                try:
                    # Align the two sequences and determine a similarity score.
                    alignment = seq_aligner.align(seq_one, seq_two)
                    score = alignment.score
                    max_len = max(len(seq_one), len(seq_two))

                    # Convert alignment score to an evolutionary distance.
                    pairwise_distance = max_len - score
                    row.append(pairwise_distance)

                except Exception as e:
                    print(f"Error during sequence alignment between {seq_ids[i]} and {seq_ids[j]}: {e}")
                    return None, None

        lower_triangular_matrix.append(row)

    return seq_ids, lower_triangular_matrix


def visualize_tree(tree):
    print("Your Phylogenetic Tree (ASCII Representation):")
    Phylo.draw_ascii(tree)

    # Tree is displayed graphically, if possible.
    try:
        plt.figure(figsize=(10, 6))
        Phylo.draw(tree)
        plt.show()
    except Exception as e:
        print(f"Visualization failed: {e}")


def user_input():
    construction_options = ["UPGMA", "NJ"]
    fasta_file = input("Enter the path to your FASTA file: ").strip()
    method = input("Enter the tree construction method (UPGMA/NJ): ").strip().upper()
    if method not in construction_options:
        print("Invalid method selected. Defaulting to UPGMA.")
        method = "UPGMA"
    return fasta_file, method


if __name__ == "__main__":
    # Get user's FASTA file and preferred tree construction method.
    selected_fasta_file, selected_method = user_input()

    # Retrieve sequences from the FASTA file and compute pairwise distances.
    sequences = parse_fasta_file(selected_fasta_file)
    sequence_ids, distance_matrix = calculate_distance_matrix(sequences)

    # Initialize class to be able to create the tree.
    constructor = DistanceTreeConstructor()

    # Convert the distance matrix into a format required for tree construction.
    dm = _DistanceMatrix(names=sequence_ids, matrix=distance_matrix)

    tree = None
    if selected_method == "UPGMA":
        # UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm is used to construct the tree.
        tree = constructor.upgma(dm)
    elif selected_method == "NJ":
        # NJ (Neighbor-Joining) algorithm is used to construct the tree.
        tree = constructor.nj(dm)

    # The resulting tree is displayed in Newick format as a textual representation.
    print("Phylogenetic Tree (Newick Format):")
    print(tree.format("newick"))

    # Visualize the tree in ASCII format and optionally as a graphical representation.
    visualize_tree(tree)










