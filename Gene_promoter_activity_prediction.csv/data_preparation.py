import pandas as pd
import os

def convert_tsv_to_csv(tsv_file="data/gene_expression.A549.biorep1.tsv", csv_file="data/expression_data.csv"):
    """
    Converts a TSV file to a CSV file and cleans up the gene IDs.

    Args:
        tsv_file (str): The path to the input TSV file.
        csv_file (str): The path where the output CSV file will be saved.

    Returns:
        pandas.DataFrame: The content of the CSV file as a DataFrame.
    """
    # Read the tab-separated file into a pandas DataFrame.
    df = pd.read_csv(tsv_file, sep="\t")
    
    # Clean up the 'gene_id' column by removing version numbers (e.g., '.1', '.2').
    df['gene_id'] = df['gene_id'].apply(lambda x: x.split('.')[0])
    
    # Announce the conversion from TSV to CSV.
    print("✅ Converting from tsv to csv for efficient usage.")
    
    # Save the DataFrame to a CSV file without including the DataFrame index.
    df.to_csv(csv_file, index=False)
    
    return df

def merge_expression_promoters(expr_file="data/expression_data.csv", prom_file="data/promoter_sequences.csv", output_file="data/merged_output.csv"):
    """
    Merges gene expression data with promoter sequences.

    Args:
        expr_file (str): The path to the expression data CSV file.
        prom_file (str): The path to the promoter sequences CSV file.
        output_file (str): The path where the merged CSV file will be saved.

    Returns:
        pandas.DataFrame: The merged data as a DataFrame.
    """
    # Read the expression and promoter data into separate DataFrames.
    df_expr = pd.read_csv(expr_file)
    df_prom = pd.read_csv(prom_file)
    
    # Merge the two DataFrames based on the 'gene_id' column.
    merged_df = pd.merge(df_expr, df_prom, on='gene_id', how='inner')
    
    # Save the merged DataFrame to a new CSV file.
    merged_df.to_csv(output_file, index=False)
    
    # Announce the successful merging of files.
    print("✅ Merging files for processing")
    
    return merged_df

def create_final_dataset(motif_file="data/promoter_motif_features.csv", onehot_file="data/promoter_onehot_100bp.csv", output_file="data/Final_dataset.csv"):
    """
    Combines motif features and one-hot encoded data to create the final dataset.

    Args:
        motif_file (str): The path to the promoter motif features CSV file.
        onehot_file (str): The path to the one-hot encoded promoter data CSV file.
        output_file (str): The path where the final dataset will be saved.

    Returns:
        pandas.DataFrame or None: The final dataset as a DataFrame, or None if input files are not found.
    """
    # Verify that both input files exist before proceeding.
    if not os.path.exists(motif_file):
        print(f"❌ Error: Input file '{motif_file}' not found.")
        return None
    if not os.path.exists(onehot_file):
        print(f"❌ Error: Input file '{onehot_file}' not found.")
        return None

    # Load the motif and one-hot encoded data into DataFrames.
    df_motif = pd.read_csv(motif_file)
    df_onehot = pd.read_csv(onehot_file)

    # Merge the two DataFrames on the 'gene_id' column.
    df = pd.merge(df_onehot, df_motif, on='gene_id', how='inner')
    
    # If 'length' and 'effective_length' columns are missing, calculate them from the promoter sequence.
    if 'length' not in df.columns and 'promoter_sequence' in df.columns:
        df['length'] = df['promoter_sequence'].apply(lambda x: len(str(x)) if pd.notna(x) else 0)
    if 'effective_length' not in df.columns and 'promoter_sequence' in df.columns:
        df['effective_length'] = df['promoter_sequence'].apply(
            lambda x: len([c for c in str(x).upper() if c in 'ACGT']) if pd.notna(x) else 0
        )

    # Remove the 'transcript_id(s)' and 'promoter_sequence' columns to avoid redundancy.
    df = df.drop(columns=["transcript_id(s)"], axis=1)
    if 'promoter_sequence' in df.columns:
        df = df.drop(columns=['promoter_sequence'])

    # Reorder columns to move 'TPM' and 'FPKM' to the end.
    columns = [col for col in df.columns if col not in ['TPM', 'FPKM']]
    if 'TPM' in df.columns:
        columns.append('TPM')
    if 'FPKM' in df.columns:
        columns.append('FPKM')
    df = df[columns]

    # Save the final dataset to a CSV file.
    df.to_csv(output_file, index=False)
    print(f"✅ Final dataset saved to '{output_file}'")
    
    return df 