def calculate_gc_content(seq):
    """Calculate the GC content of a sequence as a percentage."""
    if not seq:
        return 0.0
    seq = seq.upper()
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in seq:
     if base in counts:
        counts[base] += 1
    print("Counts:", counts)
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100

