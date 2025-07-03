import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2


def alignment(seq_1, seq_2, match=2, mismatch=-2, gap_open=-2, gap_extended=-2):
    align = pairwise2.align.globalms(seq_1, seq_2, match, mismatch, gap_open, gap_extended)
    best_alignment = align[0]
    
    print("Aligned Sequence score:", best_alignment.score,"\n")
    print("Aligned Sequence 1 :", best_alignment.seqA[:30],"........")
    print("Aligned Sequence 2 :", best_alignment.seqB[:30],"........\n")
    print("start of Aligned Sequence:", best_alignment.start)
    print("end of Aligned Sequence:", best_alignment.end)
    return best_alignment


def similarity(align):
    seq_1 = align.seqA
    seq_2 = align.seqB
    start = align.start
    end = align.end
    aligned_1 = seq_1[start:end]
    aligned_2 = seq_2[start:end]
    matches = 0
    for i in range(len(aligned_1)):
        if aligned_1[i] == aligned_2[i] and aligned_1[i] != '-':
            matches += 1
    length = end - start
    similarity = (matches / length) * 100 if length > 0 else 0
    print("The Similarity of aligned sequence is:", round(similarity,2), "%")
    return similarity


def gap_freq(align):
    seq_1 = align.seqA
    seq_2 = align.seqB
    gasp_1 = seq_1.count("-")
    gasp_2 = seq_2.count("-")
    f_gap_1 = (gasp_1 / len(seq_1)) * 100
    f_gap_2 = (gasp_2 / len(seq_2)) * 100
    print("Gap in seq_1:", round(f_gap_1,2),"%")
    print("Gap in seq_2:", round(f_gap_2,2),"%")


def conserved_reg(align):
    seq_1 = align.seqA
    seq_2 = align.seqB
    conserve_reg = []
    track_reg = []
    nucleotide_count = 0
    seq_1_index = -1  
    seq_2_index = -1  

    for i in range(len(seq_1)):
        if seq_1[i] != '-':
            seq_1_index += 1
        if seq_2[i] != '-':
            seq_2_index += 1

        
        if seq_1[i] == seq_2[i] and seq_1[i] != '-' and seq_2[i] != '-':
            if not track_reg:
                track_reg = [seq_1_index, seq_1_index, seq_2_index, seq_2_index]  
            track_reg[1] = seq_1_index 
            track_reg[3] = seq_2_index 
            nucleotide_count += 1
        else:
            
            if track_reg and nucleotide_count >= 20:
                conserve_reg.append((track_reg[0], track_reg[1], track_reg[2], track_reg[3], nucleotide_count))
            track_reg = []
            nucleotide_count = 0


    if track_reg and nucleotide_count >= 20:
        conserve_reg.append((track_reg[0], track_reg[1], track_reg[2], track_reg[3], nucleotide_count))

    
    print("Conserved regions (with 20 0r more nucleotides):")
    if conserve_reg:
        for start1, end1, start2, end2, length in conserve_reg:
            print(f"Region: Seq_1, ({start1} - {end1}) Seq_2 ({start2} - {end2}) length: {length} nucleotides")
    else:
        print("No conserved regions of 20 or more nucleotides found.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python Project.py sequence1 sequence2")
    
    file_1 = sys.argv[1]
    file_2 = sys.argv[2]
    seq2 = SeqIO.read(file_1,"fasta")
    seq1 = SeqIO.read(file_2,"fasta")
    seq_1 = seq1.seq
    seq_2 = seq2.seq
    print("Length of sequence_1:", len(seq_1))
    print("Length of sequence_2:", len(seq_2))

    align = alignment(seq_1, seq_2)
    similarity(align)
    gap_freq(align)
    conserved_reg(align)
