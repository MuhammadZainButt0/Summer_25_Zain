import pandas as pd
import requests
from tqdm import tqdm
import time

def fetch_gene_info(gene_id):
    """
    Fetches gene information, including coordinates, from the Ensembl database.

    Args:
        gene_id (str): The Ensembl ID of the gene.

    Returns:
        dict or None: A dictionary with gene information, or None if the request fails.
    """
    # The base URL for the Ensembl REST API.
    server = "https://rest.ensembl.org"
    # The specific API endpoint for looking up a gene ID.
    ext = f"/lookup/id/{gene_id}?expand=1"
    # Set the content type for the request.
    headers = {"Content-Type": "application/json"}
    
    try:
        # Make the GET request to the Ensembl API.
        r = requests.get(server + ext, headers=headers, timeout=10)
        # If the request was successful, return the JSON response.
        if r.ok:
            return r.json()
        else:
            print(f"⚠️ Failed to fetch info for {gene_id}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"⚠️ Error for {gene_id}: {e}")
        return None

def fetch_promoter_sequence(chrom, start, end, strand, species="human"):
    """
    Fetches a DNA sequence for a specific genomic region from the Ensembl database.

    Args:
        chrom (str): The chromosome number.
        start (int): The start position of the sequence.
        end (int): The end position of the sequence.
        strand (int): The strand of the DNA (1 for forward, -1 for reverse).
        species (str): The species to fetch the sequence from.

    Returns:
        str or None: The DNA sequence as a string, or None if the request fails.
    """
    server = "https://rest.ensembl.org"
    # Format the genomic region string.
    region = f"{chrom}:{start}..{end}:{strand}"
    # The API endpoint for fetching a sequence from a region.
    ext = f"/sequence/region/{species}/{region}?"
    headers = {"Content-Type": "text/plain"}
    
    try:
        r = requests.get(server + ext, headers=headers, timeout=10)
        # If successful, return the sequence text.
        if r.ok:
            return r.text.strip()
        else:
            print(f"⚠️ Failed to fetch sequence for {region}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"⚠️ Error for {region}: {e}")
        return None

def get_promoter_sequence(gene_id, upstream=1000):
    """
    Gets the promoter sequence for a given gene.

    Args:
        gene_id (str): The Ensembl ID of the gene.
        upstream (int): The number of base pairs upstream of the gene to fetch.

    Returns:
        str or None: The promoter sequence, or None if it cannot be fetched.
    """
    # First, fetch the gene's information.
    info = fetch_gene_info(gene_id)
    if info is None or 'seq_region_name' not in info:
        return None
        
    # Extract the necessary information to define the promoter region.
    chrom = info['seq_region_name']
    strand = info['strand']
    gene_start = info['start']
    gene_end = info['end']
    
    # Calculate the start and end of the promoter based on the gene's strand.
    if strand == 1:  # Forward strand
        start = max(1, gene_start - upstream)
        end = gene_start - 1
    else:  # Reverse strand
        start = gene_end + 1
        end = gene_end + upstream
        
    # Fetch the promoter sequence using the calculated coordinates.
    seq = fetch_promoter_sequence(chrom, start, end, strand)
    return seq

def fetch_promoter_sequences(gene_ids, output_file="data/promoter_sequences.csv"):
    """
    Fetches promoter sequences for a list of gene IDs and saves them to a CSV file.

    Args:
        gene_ids (list): A list of Ensembl gene IDs.
        output_file (str): The path to save the output CSV file.

    Returns:
        pandas.DataFrame: A DataFrame containing the gene IDs and their promoter sequences.
    """
    promoters = []
    # Loop through the list of gene IDs with a progress bar.
    for gid in tqdm(gene_ids, desc="Fetching promoter sequences"):
        seq = get_promoter_sequence(gid, upstream=1000)
        promoters.append({"gene_id": gid, "promoter_sequence": seq if seq else ""})
        time.sleep(0.2)  # Pause to avoid overwhelming the API server.
        
    # Convert the list of promoter data into a DataFrame and save it to a CSV.
    df = pd.DataFrame(promoters)
    df.to_csv(output_file, index=False)
    print(f"✅ Saved promoter sequences to '{output_file}'")
    
    return df 