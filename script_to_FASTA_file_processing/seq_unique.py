def get_unique_nucleotides(fasta_dict):
    """Identify unique nucleotides across all sequences using a set."""
    unique_nucleotides = set()
    for seq in fasta_dict.values():
        unique_nucleotides.update(seq.upper())
    return sorted(unique_nucleotides)