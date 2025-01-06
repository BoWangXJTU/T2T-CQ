#!/usr/bin/env python3
import argparse
from collections import defaultdict
import numpy as np
from multiprocessing import Pool, cpu_count

"""
Optimized Genome Statistics Calculation:
- Genome size, scaffold, and contig lengths.
- Base frequencies (N, G, C, T, A).
- N50/N90 values.
- Parallel processing support.

USAGE: ./genome_statistics_optimized.py -i genome.fasta -o stats.txt [-t threads]
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
                sequences[current_seq] = []
            else:
                sequences[current_seq].append(line.upper())
    # Convert lists to strings
    sequences = {k: ''.join(v) for k, v in sequences.items()}
    return sequences

def process_sequence(seq):
    """Process a single sequence to calculate length, base counts, and contig lengths."""
    length = len(seq)
    N_num = seq.count('N')
    G_num = seq.count('G')
    C_num = seq.count('C')
    T_num = seq.count('T')
    A_num = seq.count('A')
    contigs = [contig for contig in seq.split('N' * 10) if contig]
    contig_lengths = [len(contig) for contig in contigs]
    return length, N_num, G_num, C_num, T_num, A_num, contig_lengths

def calculate_statistics(sequences, threads):
    """Calculate genome statistics from sequences using parallel processing."""
    genome_size = 0
    scaffold_lengths = []
    contig_lengths = []
    N_num = G_num = C_num = T_num = A_num = 0

    # Use multiprocessing to process sequences in parallel
    with Pool(threads) as pool:
        results = pool.map(process_sequence, sequences.values())

    # Aggregate results from parallel processing
    for result in results:
        length, n, g, c, t, a, contigs = result
        genome_size += length
        scaffold_lengths.append(length)
        N_num += n
        G_num += g
        C_num += c
        T_num += t
        A_num += a
        contig_lengths.extend(contigs)

    # Sort lengths in descending order
    scaffold_lengths = np.array(scaffold_lengths)
    scaffold_lengths.sort()
    scaffold_lengths = scaffold_lengths[::-1]  # Reverse to get descending order

    contig_lengths = np.array(contig_lengths)
    contig_lengths.sort()
    contig_lengths = contig_lengths[::-1]  # Reverse to get descending order

    def calculate_Nx(lengths, x, total_size):
        """Calculate Nx value (e.g., N50, N90)."""
        target = total_size * (x / 100)
        cumulative = np.cumsum(lengths)
        idx = np.where(cumulative >= target)[0][0]
        return lengths[idx]

    # Calculate statistics
    stats = {
        "genome_size": genome_size,
        "scaffold_count": len(scaffold_lengths),
        "longest_scaffold": scaffold_lengths[0] if len(scaffold_lengths) > 0 else 0,
        "shortest_scaffold": scaffold_lengths[-1] if len(scaffold_lengths) > 0 else 0,
        "rate_of_N": N_num / genome_size if genome_size else 0,
        "rate_of_GC": (G_num + C_num) / (G_num + C_num + T_num + A_num) if (G_num + C_num + T_num + A_num) else 0,
        "scaffold_N50": calculate_Nx(scaffold_lengths, 50, genome_size),
        "scaffold_N90": calculate_Nx(scaffold_lengths, 90, genome_size),
        "contig_N50": calculate_Nx(contig_lengths, 50, genome_size),
        "contig_N90": calculate_Nx(contig_lengths, 90, genome_size),
        "sequences_1kb": np.sum(scaffold_lengths >= 1000),
        "total_length_1kb": np.sum(scaffold_lengths[scaffold_lengths >= 1000]),
        "sequences_2kb": np.sum(scaffold_lengths >= 2000),
        "total_length_2kb": np.sum(scaffold_lengths[scaffold_lengths >= 2000]),
        "sequences_3kb": np.sum(scaffold_lengths >= 3000),
        "total_length_3kb": np.sum(scaffold_lengths[scaffold_lengths >= 3000]),
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
    parser.add_argument('-t', '--threads', type=int, default=cpu_count(), help='Number of threads to use (default: all cores)')
    args = parser.parse_args()

    # Parse the input FASTA file
    sequences = parse_fasta(args.input)

    # Calculate statistics using parallel processing
    stats = calculate_statistics(sequences, args.threads)

    # Write the results to the output file
    write_output(stats, args.output)
    print(f"Genome statistics written to {args.output}")

if __name__ == '__main__':
    main()