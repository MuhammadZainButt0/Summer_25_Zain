import pandas as pd
import csv

def one_hot_encode_promoter(input_file="data/promoter_motif_features.csv", output_file="data/promoter_onehot_100bp.csv", fixed_length=100):
    """
    Converts the first 100 bases of promoter sequences into a one-hot encoded format.

    Args:
        input_file (str): The path to the input CSV file containing promoter sequences.
        output_file (str): The path to save the one-hot encoded data.
        fixed_length (int): The length of the sequence to encode.
    """
    # A dictionary to map each DNA base to its one-hot encoding.
    one_hot = {'A': ['1', '0', '0', '0'], 'C': ['0', '1', '0', '0'], 'G': ['0', '0', '1', '0'], 'T': ['0', '0', '0', '1'], 'N': ['0', '0', '0', '0']}
    
    # Read the input data.
    df = pd.read_csv(input_file)
    
    # Open the output file for writing.
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Create the header row for the CSV file.
        header = ['gene_id']
        for i in range(fixed_length):
            header.extend([f'pos{i+1}_A', f'pos{i+1}_C', f'pos{i+1}_G', f'pos{i+1}_T'])
        writer.writerow(header)
        
        # Iterate over each row in the input DataFrame.
        for _, row in df.iterrows():
            gene_id = row['gene_id']
            # Get the first 100 bases of the promoter sequence and pad with 'N' if shorter.
            seq = str(row['promoter_sequence']).upper()[:fixed_length].ljust(fixed_length, 'N')
            
            row_data = [gene_id]
            # Convert each base in the sequence to its one-hot encoding.
            for base in seq:
                row_data.extend(one_hot.get(base, one_hot['N']))
            
            # Write the one-hot encoded data to the CSV file.
            writer.writerow(row_data)
            
    print(f"✅ One-hot encoding of Promoter Sequences with fixed length {fixed_length} bases")
    
    # Return the one-hot encoded data as a DataFrame.
    return pd.read_csv(output_file) 