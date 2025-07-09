import numpy as np
import pandas as pd
import sys
import Bio
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def protein_features_analsis(fasta_file):
    try:
        data= []

        for record in SeqIO.parse(fasta_file,"fasta"):
            # sequence = record.seq
            sequence = str(record.seq).replace("-","").replace("*","")
            protein = ProteinAnalysis(str(sequence))
            features = {
                "Protein ID" : record.id,
                "Molecular_weight" : protein.molecular_weight(),
                "Isoelectric_point": protein.isoelectric_point(),
                "Aromaticity" : protein.aromaticity(),
                "Instability_index" : protein.instability_index(),
                "Gravy": protein.gravy()
            }
            data.append (features)
        
        df = pd.DataFrame(data)
        # logging.info(df.to_string())
        return df
    except:
        # logging.error(f"File not found")
        print("error in file")
        sys.exit(1)



def synthetic_label(df):
        try:
             
             labels = ["cytoplasm","mitochondria","nucleus"]
             df["labels"] = np.random.choice(labels , size=len(df))

             logging.info(f"Given_Labelled_data:\n{df}")
             return
             
        except:
              print("error in  lable file")
              sys.exit(1)