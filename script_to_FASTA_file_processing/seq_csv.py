import csv
from seq_gc import calculate_gc_content
from seq_validation import validate_seq

def write_dict_to_csv(fasta_dict, output_file):
    """Write sequence IDs, lengths, GC content, and validation results to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence_ID', 'Sequence_Length', 'GC_Content (%)', 'Validation_Status', 'Validation_Message','Sequence'])
        for seq_id, seq in fasta_dict.items():
            seq_length = len(seq)
            gc_content = calculate_gc_content(seq)
            is_valid, message = validate_seq(seq)
            writer.writerow([seq_id, seq_length, f"{gc_content:.2f}", 'Valid' if is_valid else 'Invalid', message,seq])