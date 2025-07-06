import os
import pandas as pd
from collections import defaultdict
from intervaltree import IntervalTree

def parse_gff(gff_file):
    gene_exons = defaultdict(list)
    transcript_to_gene = {}  # Map transcript_id to gene_id
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found at: {gff_file}")
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            seqid, _, feature, start, end, _, strand, _, attributes = fields
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                parts = attr.split('=', 1)
                if len(parts) == 2:
                    key, val = parts
                    attr_dict[key] = val.strip('"')
            
            if feature == 'gene':
                gene_id = attr_dict.get('gene_id')
                if gene_id:
                    continue
            elif feature == 'mRNA':
                transcript_id = attr_dict.get('transcript_id')
                parent_gene = attr_dict.get('Parent')
                if transcript_id and parent_gene and parent_gene.startswith('gene:'):
                    transcript_to_gene[transcript_id] = parent_gene.replace('gene:', '')
                    continue
            elif feature != 'exon':
                continue
            
            parent = attr_dict.get('Parent')
            if parent and parent.startswith('transcript:'):
                transcript_id = parent.replace('transcript:', '')
                gene_id = transcript_to_gene.get(transcript_id)
                if gene_id:
                    gene_exons[gene_id].append({
                        'seqid': seqid,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand
                    })
    return gene_exons

def parse_snps(vcf_file):
    snps = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])
            snps.append({'chrom': chrom, 'pos': pos})
        print(f"Parsed {len(snps)} SNPs")
    return snps

def map_snps_to_exons(snps, gene_exons):
    chrom_trees = defaultdict(IntervalTree)
    for gene, exons in gene_exons.items():
        for exon in exons:
            chrom = exon['seqid']
            start = exon['start']
            end = exon['end'] + 1
            chrom_trees[chrom][start:end] = {'gene': gene, 'exon': exon}
    
    snp_to_exons = defaultdict(list)
    for snp in snps:
        chrom = snp['chrom']
        pos = snp['pos']
        if chrom in chrom_trees:
            overlapping_intervals = chrom_trees[chrom][pos]
            for interval in overlapping_intervals:
                snp_to_exons[(chrom, pos)].append(interval.data)
    print(f"Mapped {len(snp_to_exons)} SNPs to exons")
    return snp_to_exons

def calculate_snp_density(gene_exons, snp_mappings):
    gene_stats = defaultdict(lambda: {'length': 0, 'snp_count': 0})
    for gene, exons in gene_exons.items():
        for exon in exons:
            length = exon['end'] - exon['start'] + 1
            gene_stats[gene]['length'] += length
    for (chrom, pos), mappings in snp_mappings.items():
        for mapping in mappings:
            gene = mapping['gene']
            gene_stats[gene]['snp_count'] += 1
    density_results = {}
    for gene, stats in gene_stats.items():
        length_kb = stats['length'] / 1000.0
        snp_count = stats['snp_count']
        density = snp_count / length_kb if length_kb > 0 else 0
        density_results[gene] = {
            'total_exonic_length': stats['length'],
            'snp_count': snp_count,
            'density': density
        }
    return density_results

# Example usage
gff_file = "Saccharomyces_cerevisiae.R64-1-1.61.gff3"  # Updated to match likely file name
vcf_file = "saccharomyces_cerevisiae.vcf"  # Verify this path
gene_exons = parse_gff(gff_file)
snps = parse_snps(vcf_file)
snp_mappings = map_snps_to_exons(snps, gene_exons)
density_results = calculate_snp_density(gene_exons, snp_mappings)

# Output results
with open("SNP_density.csv", "w") as f:
    f.write("Gene,Total_Exonic_Length,SNP_Count,Density_SNPs_per_kb\n")
    for gene, stats in density_results.items():
        line = f"{gene},{stats['total_exonic_length']},{stats['snp_count']},{stats['density']:.2f}\n"
        f.write(line)
    print("Total Exonic Regions:", len(density_results))
    print("Results are saved to P_SNP_density.csv file")