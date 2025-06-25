def validate_seq(seq):
    """Validate if a sequence contains only valid DNA nucleotides (A, T, C, G)."""
    seq = seq.upper()
    nucleotides = {"A", "T", "C", "G"}
    invalid_bases = [base for base in seq if base not in nucleotides]
    if invalid_bases:
        return False, f"Invalid nucleotides found: {invalid_bases}"
    return True, "Valid DNA sequence"