import pandas as pd
import numpy as np
from collections import Counter
from itertools import product
import re
import math

# The following functions calculate various features from a DNA sequence.

def gc_content(seq):
    """Calculates the percentage of G and C bases in a sequence."""
    if pd.isna(seq) or len(seq) == 0:
        return 0.0
    seq = str(seq).upper()
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 2) if len(seq) > 0 else 0.0

def at_content(seq):
    """Calculates the percentage of A and T bases in a sequence."""
    if pd.isna(seq) or len(seq) == 0:
        return 0.0
    seq = seq.upper()
    a = seq.count('A')
    t = seq.count('T')
    return round(100 * (a + t) / len(seq), 2) if len(seq) > 0 else 0.0

def cpg_density(seq):
    """Calculates the density of CpG dinucleotides per 100 base pairs."""
    if pd.isna(seq) or len(seq) == 0:
        return 0.0
    seq = seq.upper()
    cpgs = seq.count("CG")
    return round((cpgs / len(seq)) * 100, 2) if len(seq) > 0 else 0.0

def shannon_entropy(seq):
    """Calculates the Shannon entropy, a measure of sequence complexity."""
    if pd.isna(seq) or len(seq) == 0:
        return 0.0
    seq = seq.upper()
    bases = ['A', 'C', 'G', 'T']
    # Calculate the frequency of each base.
    freq = [seq.count(base) / len(seq) for base in bases]
    # Calculate the entropy.
    entropy = -sum(p * math.log2(p) for p in freq if p > 0)
    return round(entropy, 4)

def kmer_frequencies(seq, kmer_len=3):
    """Calculates the frequency of all possible k-mers (short DNA words)."""
    if pd.isna(seq) or len(seq) < kmer_len:
        # Return a dictionary of zeros if the sequence is too short.
        return {f'kmer_{"".join(p)}': 0 for p in product('ACGT', repeat=kmer_len)}
    
    seq = str(seq).upper()
    total_kmers = len(seq) - kmer_len + 1
    # Get all k-mers from the sequence.
    kmers = [seq[i:i+kmer_len] for i in range(total_kmers)]
    # Count the occurrences of each k-mer.
    kmer_counts = Counter(kmers)
    
    # Normalize the counts to get frequencies.
    for k in kmer_counts:
        kmer_counts[k] /= total_kmers
        
    all_kmers = [''.join(p) for p in product('ACGT', repeat=kmer_len)]
    # Return a dictionary with frequencies for all possible k-mers.
    return {f'kmer_{kmer}': kmer_counts.get(kmer, 0) for kmer in all_kmers}

def dinucleotide_freq(seq):
    """Calculates the frequency of all possible dinucleotides (pairs of bases)."""
    dinucleotides = [''.join(p) for p in product('ACGT', repeat=2)]
    if pd.isna(seq) or len(seq) < 2:
        return {f"freq_{d}": 0.0 for d in dinucleotides}
        
    seq = seq.upper()
    counts = Counter([seq[i:i+2] for i in range(len(seq)-1)])
    total = sum(counts.values())
    
    return {f"freq_{d}": round(100 * counts.get(d, 0) / total, 2) if total > 0 else 0.0 for d in dinucleotides}

def search_motif(seq, motif, motif_name):
    """Searches for a specific DNA motif in the last 100 bp of a sequence."""
    if pd.isna(seq) or len(seq) == 0:
        return pd.Series([0, -1, 0], index=[f"{motif_name}_presence", f"{motif_name}_position", f"{motif_name}_count"])
        
    seq = str(seq).upper()
    # Search only in the last 100 base pairs.
    window = seq[-100:] if len(seq) > 100 else seq
    pattern = re.compile(motif)
    matches = list(pattern.finditer(window))
    
    presence = int(len(matches) > 0)
    position = matches[0].start() if matches else -1
    count = len(matches)
    
    return pd.Series([presence, position, count], index=[f"{motif_name}_presence", f"{motif_name}_position", f"{motif_name}_count"])

def extract_sequence_features(input_file="data/merged_output.csv", output_file="data/feature_extracted_output.csv"):
    """
    Extracts a variety of sequence-based features and adds them to the dataset.
    """
    df = pd.read_csv(input_file)
    df['promoter_sequence'] = df['promoter_sequence'].fillna('')
    
    # Apply each feature extraction function to the promoter sequence column.
    df['GC_Content'] = df['promoter_sequence'].apply(gc_content)
    df['AT_content_%'] = df['promoter_sequence'].apply(at_content)
    df['CpG_density_per_100bp'] = df['promoter_sequence'].apply(cpg_density)
    df['promoter_length'] = df['promoter_sequence'].apply(len)
    df['shannon_entropy'] = df['promoter_sequence'].apply(shannon_entropy)
    
    # Get k-mer and dinucleotide frequencies and add them as new columns.
    kmer_df = pd.DataFrame(df['promoter_sequence'].apply(kmer_frequencies, kmer_len=3).tolist())
    dinuc_df = pd.DataFrame(df['promoter_sequence'].apply(dinucleotide_freq).apply(pd.Series))
    df = pd.concat([df, kmer_df, dinuc_df], axis=1)
    
    df.to_csv(output_file, index=False)
    print("✅ Sequence features extracted Successfully'")
    return df

def extract_motif_features(input_file="data/feature_extracted_output.csv", output_file="data/promoter_motif_features.csv"):
    """
    Extracts features based on the presence and location of common regulatory motifs.
    """
    df = pd.read_csv(input_file)
    df['promoter_sequence'] = df['promoter_sequence'].fillna('')
    
    # Special handling for the TATA box motif.
    def extract_tata_features(seq):
        # Search for the TATA box in a specific window upstream of the gene.
        window_start = max(len(seq) - 50, 0)
        window_end = max(len(seq) - 20, 0)
        search_window = seq[window_start:window_end]
        motif = "TATAAA"
        
        presence = 1 if motif in search_window else 0
        position = search_window.find(motif) if presence else -1
        count = search_window.count(motif)
        
        return pd.Series([presence, position, count], index=["tata_presence", "tata_position", "tata_count"])
        
    df[["tata_presence", "tata_position", "tata_count"]] = df["promoter_sequence"].apply(extract_tata_features)
    
    # A dictionary of other motifs to search for.
    motifs = {
        "caat": "CCAAT",
        "gc": "GGGCGG",
        "inr": r"[CT]{2}A[AT][CT]{2}"  # Initiator element (uses regular expression)
    }
    
    # Search for each motif and add the results as new columns.
    for name, pattern in motifs.items():
        motif_features = df["promoter_sequence"].apply(search_motif, args=(pattern, name))
        df = pd.concat([df, motif_features], axis=1)
        
    df.to_csv(output_file, index=False)
    print("✅ Motif features extracted Successfully............................")
    return df 