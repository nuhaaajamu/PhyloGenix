from Bio import SeqIO


def parse_fasta_file(file_path):
    parsed_sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        # Retrieve the sequence identifier and corresponding sequence
        parsed_sequences[record.id] = str(record.seq)
    return parsed_sequences
