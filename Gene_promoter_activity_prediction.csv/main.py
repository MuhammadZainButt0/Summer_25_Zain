import os
from data_preparation import convert_tsv_to_csv, merge_expression_promoters, create_final_dataset
from sequence_fetcher import fetch_gene_info,fetch_promoter_sequence,get_promoter_sequence,fetch_promoter_sequences
from feature_extraction import extract_sequence_features, extract_motif_features
from one_hot_encoding import one_hot_encode_promoter

from model_training import train_evaluate_model, delete_csvs_except_final

def main():
    """
    Main function to run the entire data processing and model training pipeline.
    """
    # A list of Ensembl gene IDs to be used in the analysis.
    gene_ids = [
        "ENSG00000000003", "ENSG00000000005", "ENSG00000000419", "ENSG00000000457", "ENSG00000000460",
        "ENSG00000000938", "ENSG00000000971", "ENSG00000001036", "ENSG00000001084", "ENSG00000001167",
        "ENSG00000001460", "ENSG00000001461", "ENSG00000001497", "ENSG00000001561", "ENSG00000001617",
        "ENSG00000001626", "ENSG00000001629", "ENSG00000001630", "ENSG00000001631", "ENSG00000002016",
        "ENSG00000002079", "ENSG00000002330", "ENSG00000002549", "ENSG00000002586", "ENSG00000002587",
        "ENSG00000002726", "ENSG00000002745", "ENSG00000002746", "ENSG00000002822", "ENSG00000002834",
        "ENSG00000002919", "ENSG00000002933", "ENSG00000003056", "ENSG00000003096", "ENSG00000003137",
        "ENSG00000003147", "ENSG00000003249", "ENSG00000003393", "ENSG00000003400", "ENSG00000003402",
        "ENSG00000003436", "ENSG00000003509", "ENSG00000003756", "ENSG00000003987", "ENSG00000003989",
        "ENSG00000004059", "ENSG00000004139", "ENSG00000004142", "ENSG00000004399", "ENSG00000004455",
        "ENSG00000004468", "ENSG00000004478", "ENSG00000004487", "ENSG00000004534", "ENSG00000004660",
        "ENSG00000004700", "ENSG00000004766", "ENSG00000004776", "ENSG00000004777", "ENSG00000004779",
        "ENSG00000004799", "ENSG00000004809", "ENSG00000004838", "ENSG00000004846", "ENSG00000004848",
        "ENSG00000004864", "ENSG00000004866", "ENSG00000004897", "ENSG00000004939", "ENSG00000004948",
        "ENSG00000004961", "ENSG00000004975", "ENSG00000005001", "ENSG00000005007", "ENSG00000005020",
        "ENSG00000005022", "ENSG00000005059", "ENSG00000005073", "ENSG00000005075", "ENSG00000005100",
        "ENSG00000005102", "ENSG00000005108", "ENSG00000005156", "ENSG00000005175", "ENSG00000005187",
        "ENSG00000005189", "ENSG00000005194", "ENSG00000005206", "ENSG00000005238", "ENSG00000005243",
        "ENSG00000005249", "ENSG00000005302", "ENSG00000005339", "ENSG00000005379", "ENSG00000005381",
        "ENSG00000005421", "ENSG00000005436", "ENSG00000005448", "ENSG00000005469", "ENSG00000005471"
    ]

    # Path to the initial gene expression data file.
    tsv_file = "data/gene_expression.A549.biorep1.tsv"
    
    # Check if the input TSV file exists before running the pipeline.
    if not os.path.exists(tsv_file):
        print(f"❌ Error: Input file '{tsv_file}' not found. Please ensure it exists in the 'data' directory.")
        return

    # The following steps execute the entire pipeline in sequence.
    
    # 1. Convert the raw TSV data to a more usable CSV format.
    convert_tsv_to_csv()
    
    # 2. Fetch promoter sequences for the specified genes. 
    
    fetch_promoter_sequences(gene_ids)
    
    # 3. Merge the expression data with the promoter sequences.
    merge_expression_promoters()
    
    # 4. Extract basic sequence features (e.g., GC content).
    extract_sequence_features()
    
    # 5. Extract features related to transcription factor binding motifs.
    extract_motif_features()
    
    # 6. Convert promoter sequences into a numerical format (one-hot encoding).
    one_hot_encode_promoter()
    
    # 7. Create the final, consolidated dataset for model training.
    create_final_dataset()
    
    # 8. Train a machine learning model and evaluate its performance.
    train_evaluate_model()
    
    # 9. Clean up intermediate CSV files, keeping only the final results.
    delete_csvs_except_final()

# This standard Python construct ensures that the main() function is called only when the script is executed directly.
if __name__ == "__main__":
    main() 
