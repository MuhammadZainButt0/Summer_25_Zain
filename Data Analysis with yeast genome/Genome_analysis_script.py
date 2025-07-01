import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import EcoRI
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px





def valid(seq_seq):
    try:
        for i in seq_seq:
            nucleotides = ["A","T","C","G"]
            if i not in nucleotides:
                print("\n\nSequence ID:",seq_id)
                print("Status: Invalid Sequence")
                break
            else:
                print("\n\nSequence ID:",seq_id)
                print("Status: Valid Sequence")
                break
    except:
     print("Validation is not processed")        
    return seq_seq






def GC_content(seq_seq):
  try:
    GC_content = gc_fraction(seq_seq)*100
    return GC_content
  except:
     print("No, GC content detected.")
     return None
  




def codon_freq(seq_seq):
   try:
        codons = []
        for i in range(0, len(seq_seq) - 2, 3): 
            codon = str(seq_seq[i:i+3])
            if len(codon) == 3:
                codons.append(codon)
        codon_counts = {}
        for codon in codons:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
        total_codons = sum(codon_counts.values())
        codon_freqs = {codon: count/total_codons for codon, count in codon_counts.items()}
        return codon_freqs
   except:
      print("Sequence Length is zero.")
      return {}







def RNA_Prot(seq_seq):
   try: 
    stop_codons = ["TAA", "TAG", "TGA"]
    for i in range(0, len(seq_seq) - 2, 3): 
        codon = str(seq_seq[i:i+3])
        if codon in stop_codons:
            print(f"Stop codon {codon} found at position {i}")
            Rna = seq_seq[:i].transcribe()    
            prot= Rna[:i].translate()
            print ("DNA sequence (until stop codon):",seq_seq[:i])
            print ("RNA sequence (until stop codon):",Rna)
            print ("Amino Acids Sequence (until stop codon):",prot)
            return str(Rna), str(prot)
    print("There is no Stop codon in sequence.")
    Rna = seq_seq.transcribe()
    prot = Rna.translate()
    return str(Rna), str(prot)
   except:
       print("There is no Stop codon in sequence.")
       return 





def res_sites(seq_seq):
   try:
    E_sites = EcoRI.search(seq_seq)
    total = len(E_sites)
    return total
   except:
      print("There is no EcoRI site in this sequence.")
      return 0
   




if __name__ == "__main__":
    seq = SeqIO.parse("Saccharomyces_cerevisiae.fa", "fasta")
    dict_data = []  
    sequence_ids = []  
    gc_contents = []  
    seq_lengths = []  
    ecoRI_sites = []  
    first_codon_freqs = None 

    for i, record in enumerate(seq):
        seq_id = record.id
        seq_seq = record.seq
        seq_seq = seq_seq.upper()

        valid(seq_seq)
        GC_con = GC_content(seq_seq)
        codon_freqs = codon_freq(seq_seq)
        rna_seq, prot_seq = RNA_Prot(seq_seq)
        endo_sides = res_sites(seq_seq)

        result = {
            "Sequence_ID": seq_id,
            "Length": len(seq_seq),
            "GC_Content": round(GC_con, 2) if GC_con is not None else None,
            "EcoRI_Sites": endo_sides,
            "RNA": rna_seq[:50],
            "Protein": prot_seq[:50],
            "Codon_Usage": codon_freqs
        }
        dict_data.append(result)

        
        if GC_con is not None:
            sequence_ids.append(seq_id)
            gc_contents.append(round(GC_con, 2))
            seq_lengths.append(len(seq_seq))
            ecoRI_sites.append(endo_sides)
        if i == 0: 
            first_codon_freqs = codon_freqs

        print("Number of codons in Sequence:", len(codon_freqs)) 
        print("GC content of Sequence is:", round(GC_con, 2))
        print("Sequence length is:", len(seq_seq))
        print("EcoRI sites in sequence are:", endo_sides)

    # Save results to CSV
    df = pd.DataFrame(dict_data)
    df.to_csv("sequence_analysis.csv", index=False)

    # 1. Matplotlib Bar Plot: GC Content
    fig, ax = plt.subplots()
    bar_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple'] * (len(sequence_ids) // 5 + 1)
    ax.bar(sequence_ids, gc_contents, label=sequence_ids, color=bar_colors[:len(sequence_ids)])
    ax.set_ylabel('GC Content (%)')
    ax.set_title('GC Content by Sequence')
    ax.legend(title='Sequence ID')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("gc_content_bar.png")
    plt.close()

    # 2. Seaborn Heatmap: Codon Usage Frequencies of  First Sequence
    if first_codon_freqs:
        codons = sorted(first_codon_freqs.keys())
        codon_matrix = [[first_codon_freqs.get(codon, 0) for codon in codons]]
        codon_df = pd.DataFrame(codon_matrix, columns=codons, index=[sequence_ids[0]])
        plt.figure(figsize=(12, 2))
        sns.heatmap(codon_df, annot=True, fmt=".3f", cmap="YlGnBu")
        plt.title(f"Codon Usage Frequencies for Sequence {sequence_ids[0]}")
        plt.ylabel("Sequence ID")
        plt.xlabel("Codon")
        plt.tight_layout()
        plt.savefig("codon_heatmap.png")
        plt.close()

    # 3. Matplotlib Histogram: Sequence Lengths
    plt.figure()
    plt.hist(seq_lengths, bins=10, edgecolor='black')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.title('Distribution of Sequence Lengths')
    plt.tight_layout()
    plt.savefig("length_histogram.png")
    plt.close()

    # 4. Plotly Scatter Plot: GC Content vs. EcoRI Sites
    scatter_df = pd.DataFrame({
        'Sequence_ID': sequence_ids,
        'GC_Content': gc_contents,
        'EcoRI_Sites': ecoRI_sites,
        'Length': seq_lengths
    })
    fig = px.scatter(
        scatter_df,
        x='GC_Content',
        y='EcoRI_Sites',
        size='Length',
        color='Length',
        text='Sequence_ID',
        title='GC Content vs. EcoRI Sites',
        labels={'GC_Content': 'GC Content (%)', 'EcoRI_Sites': 'EcoRI Sites', 'Length': 'Sequence Length'}
    )
    fig.update_traces(textposition='top center')
    fig.write_html("restriction_scatter.html")