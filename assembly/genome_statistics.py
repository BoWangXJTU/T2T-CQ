#!/usr/bin/env python3
import argparse
from collections import defaultdict

"""
Statistics Calculation:
Genome size, scaffold, and contig lengths.
Base frequencies (N, G, C, T, A).
N50/N90 values.

USAGE: ./genome_statistics.py -i genome.fasta -o stats.txt

"""

def parse_fasta(file):
    """Parse a FASTA file and return sequences as a dictionary."""
    sequences = {}
    current_seq = None
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq = line[1:].split()[0]
                sequences[current_seq] = ''
            else:
                sequences[current_seq] += line.upper()
    return sequences

def calculate_statistics(sequences):
    """Calculate genome statistics from sequences."""
    genome_size = 0
    scaffold_lengths = []
    contig_lengths = []
    N_num = G_num = C_num = T_num = A_num = 0
    
    for identifier, seq in sequences.items():
        length = len(seq)
        genome_size += length
        scaffold_lengths.append(length)
        
        # Base counts
        N_num += seq.count('N')
        G_num += seq.count('G')
        C_num += seq.count('C')
        T_num += seq.count('T')
        A_num += seq.count('A')
        
        # Contigs (split by stretches of 10 or more Ns)
        contigs = [contig for contig in seq.split('N' * 10) if contig]
        contig_lengths.extend(len(contig) for contig in contigs)
    
    scaffold_lengths.sort(reverse=True)
    contig_lengths.sort(reverse=True)
    
    def calculate_Nx(lengths, x, total_size):
        target = total_size * (x / 100)
        cumulative = 0
        for length in lengths:
            cumulative += length
            if cumulative >= target:
                return length
        return None
    
    stats = {
        "genome_size": genome_size,
        "scaffold_count": len(scaffold_lengths),
        "longest_scaffold": scaffold_lengths[0] if scaffold_lengths else 0,
        "shortest_scaffold": scaffold_lengths[-1] if scaffold_lengths else 0,
        "rate_of_N": N_num / genome_size if genome_size else 0,
        "rate_of_GC": (G_num + C_num) / (G_num + C_num + T_num + A_num) if (G_num + C_num + T_num + A_num) else 0,
        "scaffold_N50": calculate_Nx(scaffold_lengths, 50, genome_size),
        "scaffold_N90": calculate_Nx(scaffold_lengths, 90, genome_size),
        "contig_N50": calculate_Nx(contig_lengths, 50, genome_size),
        "contig_N90": calculate_Nx(contig_lengths, 90, genome_size),
        "sequences_1kb": sum(1 for l in scaffold_lengths if l >= 1000),
        "total_length_1kb": sum(l for l in scaffold_lengths if l >= 1000),
        "sequences_2kb": sum(1 for l in scaffold_lengths if l >= 2000),
        "total_length_2kb": sum(l for l in scaffold_lengths if l >= 2000),
        "sequences_3kb": sum(1 for l in scaffold_lengths if l >= 3000),
        "total_length_3kb": sum(l for l in scaffold_lengths if l >= 3000),
    }
    return stats

def write_output(stats, output_file):
    """Write statistics to an output file."""
    with open(output_file, 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")

def main():
    parser = argparse.ArgumentParser(description='Calculate genome statistics from a FASTA file.')
    parser.add_argument('-i', '--input', required=True, help='Input genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output statistics file')
    args = parser.parse_args()
    
    sequences = parse_fasta(args.input)
    stats = calculate_statistics(sequences)
    write_output(stats, args.output)
    print(f"Genome statistics written to {args.output}")

if __name__ == '__main__':
    main()
