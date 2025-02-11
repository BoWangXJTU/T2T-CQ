import argparse
from Bio import SeqIO
import os

# Define a function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Fill the region specified by -b in -i1 using the genome sequence from -i2')
    parser.add_argument('-i1', required=True, help='First genome file (FASTA format)')
    parser.add_argument('-i2', required=True, help='Second genome file (FASTA format)')
    parser.add_argument('-b', required=True, help='Region, e.g., 20-50')
    parser.add_argument('-o', required=True, help='Output file')
    parser.add_argument('-fold', type=int, default=1, help='Repeat the sequence from -i2 multiple times for filling')
    return parser.parse_args()

# Parse the region
def parse_region(region):
    try:
        start, end = map(int, region.split('-'))
        if start >= end or start < 1:
            raise ValueError
        return start, end
    except ValueError:
        raise ValueError('Invalid region format. Please use the correct format (e.g., 20-50)')

# Replace the sequence in the specified region
def replace_region(seq1, seq2, start, end):
    return seq1[:start-1] + seq2 + seq1[end:]

# Get the header for the output file
def get_output_header(output_file):
    return os.path.splitext(os.path.basename(output_file))[0]

# Repeat the sequence
def repeat_sequence(seq, fold):
    return seq * fold

# Main function
def main():
    args = parse_args()
    start, end = parse_region(args.b)
    
    # Read genome files
    seq1 = next(SeqIO.parse(args.i1, 'fasta')).seq
    seq2 = next(SeqIO.parse(args.i2, 'fasta')).seq
    
    # Repeat the sequence from -i2
    seq2 = repeat_sequence(seq2, args.fold)
    
    # Replace the sequence in the specified region
    filled_seq = replace_region(seq1, seq2, start, end)
    
    # Get the header for the output file
    header = get_output_header(args.o)
    
    # Write to a new FASTA file
    with open(args.o, 'w') as output:
        output.write(f'>{header}\n{filled_seq}\n')
    
    print(f'Filling completed. The result has been saved to {args.o}')

if __name__ == '__main__':
    main()
