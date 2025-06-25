from seq_to_dict import read_fasta_to_dict
from seq_gc import calculate_gc_content
from seq_validation import validate_seq
from seq_unique import get_unique_nucleotides
from seq_csv import write_dict_to_csv

def analyze_fasta_sequences(file_path, output_file):
    """Main function to read, analyze, and save FASTA sequences."""
    fasta_dict = read_fasta_to_dict(file_path)
    
    # Print analysis results for each sequence
    print(f"Analyzing sequences in {file_path}")
    for seq_id, seq in fasta_dict.items():
        seq_length = len(seq)
        gc_content = calculate_gc_content(seq)
        is_valid, message = validate_seq(seq)
        print(f"Sequence ID: {seq_id}")
        print(f"Sequence Length: {seq_length}")
        print(f"GC Content: {gc_content:.2f}%")
        print(f"Validation: {message}\n")
    
    # Print unique nucleotides across all sequences
    unique_nucleotides = get_unique_nucleotides(fasta_dict)
    print(f"Unique nucleotides across all sequences: {unique_nucleotides}")
    
    # Write results to CSV
    write_dict_to_csv(fasta_dict, output_file)
    print(f"Results saved to {output_file}")