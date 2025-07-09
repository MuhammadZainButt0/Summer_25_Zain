import logging
from model import model
from extract_features import protein_features_analsis,synthetic_label
logging.basicConfig(level=logging.INFO,filename="Result.log",format='%(asctime)s - %(levelname)s - %(message)s')


def main():
  try:  
    fasta_file = "proteins.fasta"
    
    """Parsing and analyzing features of proteins """
    df = protein_features_analsis(fasta_file)
    
    """Assigning labels to the features of proteins"""
    synthetic_label(df)

    """Model Implementation and Evaluation"""
    model(df)
    print("Successfully, results saved to log file.")
  except:
     "File not found."


if __name__ == "__main__":
    main()
    
