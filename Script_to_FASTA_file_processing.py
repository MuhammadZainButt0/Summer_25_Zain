import csv

def validate_seq(seq):
    """Validate if a sequence contains only valid DNA nucleotides (A, T, C, G)."""
    seq = seq.upper()
    nucleotides = {"A", "T", "C", "G"}
    for i in seq:
        if i not in nucleotides:
            return False, "DNA Sequence can't processed"
    return True, "DNA Sequence can processed"

def calculate_gc_content(seq):
    """Calculate the GC content of a sequence as a percentage."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0.0

def get_unique_nucleotides(fasta_dict):
    """Identify unique nucleotides across all sequences using a set."""
    unique_nucleotides = set()
    for seq in fasta_dict.values():
        unique_nucleotides.update(seq.upper())
    return unique_nucleotides

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
    
    return fasta_dict

def write_dict_to_csv(fasta_dict, output_file):
    """Write sequence IDs, lengths, GC content, and validation results to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence_ID', 'Sequence_Length', 'GC_Content (%)', 'Validation_Status', 'Validation_Message'])
        for ID, seq in fasta_dict.items():
            seq_length = len(seq)
            gc_content = calculate_gc_content(seq)
            is_valid, message = validate_seq(seq)
            writer.writerow([ID, seq_length, f"{gc_content:.2f}", 'Valid' if is_valid else 'Invalid', message])

def analyze_fasta_sequences(file_path, output_file):
    """Main function to read, analyze, and save FASTA sequences."""
    fasta_dict = read_fasta_to_dict(file_path)
    
    # Print analysis results for each sequence
    print("Analyzing sequences in", file_path)
    for ID, seq in fasta_dict.items():
        seq_length = len(seq)
        gc_content = calculate_gc_content(seq)
        is_valid, message = validate_seq(seq)
        print(f"Sequence ID: {ID}")
        print(f"Sequence Length: {seq_length}")
        print(f"GC Content: {gc_content:.2f}%")
        print(f"Validation: {message}\n")
    
    # Print unique nucleotides across all sequences
    unique_nucleotides = get_unique_nucleotides(fasta_dict)
    print(f"Unique nucleotides across all sequences: {sorted(unique_nucleotides)}")
    
    # Write results to CSV
    write_dict_to_csv(fasta_dict, output_file)
    print(f"Results saved to {output_file}")

# Usage with your file
if __name__ == "__main__":
    file_path = "sequence.fasta"
    output_file = "fasta_analysis.csv"
    analyze_fasta_sequences(file_path, output_file)