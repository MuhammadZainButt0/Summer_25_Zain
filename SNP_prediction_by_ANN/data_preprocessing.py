import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder

def load_and_preprocess_data(path="Update.csv"):
    df = pd.read_csv(path)
    # Encode Chrom
    df['Chrom'] = df['Chrom'].replace({'X': 23, 'Y': 24, 'MT': 25}).astype(int)
    # One-hot encode Ref and Alt
    ref_alt_encoded = pd.get_dummies(df[['Ref', 'Alt']], prefix=['Ref', 'Alt'])
    # Encode Variant_Type
    variant_encoder = LabelEncoder()
    df['Variant_Type_enc'] = variant_encoder.fit_transform(df['Variant_Type'])
    # Numeric features
    numeric_features = df[['Chrom','Shannon_Entropy', 'GC_Content(%)', 'AT_Content(%)','Variant_Type_enc']]
    # Combine all features except Gene and Pos
    features = pd.concat([numeric_features, ref_alt_encoded], axis=1)
    # Ensure all features are numeric
    features = features.apply(pd.to_numeric, errors='coerce').fillna(0)
    # One-hot encode ContextSeq (flattened)
    def one_hot_sequence(seq):
        mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
        seq = seq.upper()
        return [bit for base in seq for bit in mapping.get(base, [0, 0, 0, 0])]
    seq_encoded = df['ContextSeq'].apply(one_hot_sequence)
    seq_features = np.stack(seq_encoded.values)
    # Concatenate all features and cast to float32
    X = np.concatenate([features.values, seq_features], axis=1).astype(np.float32)
    y = df['Label'].values.astype(int)
    return X, y, features.columns.tolist(), seq_features.shape[1] 