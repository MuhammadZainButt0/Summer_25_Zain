import sys
import seq_analyze

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python main.py <input_fasta_file> <output_csv_file>")
    
    file_path = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        seq_analyze.analyze_fasta_sequences(file_path, output_file)
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{file_path}' not found!")
    except Exception as e:
        sys.exit(f"Unexpected error: {e}")