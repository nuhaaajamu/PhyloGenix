from Bio import SeqIO
from Bio.Align import PairwiseAligner


def parse_fasta_file(file_path):
    parsed_sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        # Retrieve the sequence identifier and corresponding sequence
        parsed_sequences[record.id] = str(record.seq)
    return parsed_sequences


def calculate_distance_matrix(parsed_sequences):
    seq_aligner = PairwiseAligner()
    seq_ids = list(parsed_sequences.keys())
    id_num = len(seq_ids)
    distance_matrix = []

    # Initialize a 2-D list that will hold calculated distances between pairs of sequences.
    for i in range(id_num):
        row = []
        for j in range (id_num):
            row.append(0)
        distance_matrix.append(row)

    for x in range(id_num):
        for y in range(x + 1, id_num):
            seq_one = parsed_sequences[seq_ids[x]]
            seq_two = parsed_sequences[seq_ids[y]]

            # Pairs of sequences are aligned and a number is assigned based on number of matches, mismatches, and gaps.
            alignment = seq_aligner.align(seq_one, seq_two)

            # The numerical value that was previously scored is retrieved.
            score = alignment.score

            max_len = max(len(parsed_sequences[seq_ids[x]]), len(parsed_sequences[seq_ids[y]]))

            # Determine evolutionary distance between sequences
            pairwise_distance = max_len - score

            # Store the distance symmetrically in the 2-D list.
            distance_matrix[x][y] = pairwise_distance
            distance_matrix[y][x] = pairwise_distance

    return seq_ids, distance_matrix









