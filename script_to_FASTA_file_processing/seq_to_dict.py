def read_fasta_to_dict(file_path):
    """Read FASTA file and store sequences in a dictionary with IDs as keys."""
    fasta_dict = {}
    current_id = None
    current_seq = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Check for header line
                if current_id is not None:  # Save previous sequence if exists
                    fasta_dict[current_id] = ''.join(current_seq)
                    current_seq = []
                current_id = line[1:]  # Store ID without '>'
            else:
                current_seq.append(line)  # Collect sequence lines
    
    # Save the last sequence
    if current_id is not None and current_seq:
        fasta_dict[current_id] = ''.join(current_seq)
    
    if not fasta_dict:
        raise ValueError("No valid sequences found in FASTA file")
    
    return fasta_dict