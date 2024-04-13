import argparse

from rotate_functions import *

parser = argparse.ArgumentParser(description="Shift circular FASTA sequences based on a reference or a position.")
parser.add_argument("-i", "--input", type=str, help="Input FASTA file. Use '-' for STDIN.")
parser.add_argument("-r", "--ref", type=str, help="Reference FASTA file.")
parser.add_argument("-p", "--pos", type=int, help="Position to shift the sequence.")
parser.add_argument("-d", "--max_dist", type=int, help="Maximum distance.")
# parser.add_argument("-lcs", "--longest_common_subsequence", action='store_true', help="Position to shift the sequence.")
args = parser.parse_args()

# Read input sequences
with open(args.input, 'r') as file_handle:
    fasta_sequences = fasta_to_dict(file_handle)
    
# Process sequences based on the reference or position
if args.ref:
    if args.max_dist:
        with open(args.ref, 'r') as ref_file:
            reference_sequences = fasta_to_dict(ref_file) # format e.g. {'reference': 'GATCGA'}
            reference_sequence = next(iter(reference_sequences.values())) # format e.g GATCGA
            for seq_id, sequence in fasta_sequences.items():
                try:
                    shifted_sequence = shift_sequence(sequence, ref=reference_sequence, max_dist=args.max_dist)
                    print(create_fasta(seq_id, shifted_sequence))
                except ValueError as e:
                    print(f"Error for {seq_id}: {str(e)}")
    else:
        raise ValueError("Please specify the max dist!")
elif args.pos:
    for seq_id, sequence in fasta_sequences.items():
        shifted_sequence = shift_sequence(sequence, pos=args.pos)
        print(create_fasta(seq_id, shifted_sequence))
# elif args.longest_common_subsequence:
#     lcs = find_longest_common_substring(*fasta_sequences.values())
#     for seq_id, sequence in fasta_sequences.items():
#         shifted_sequence = shift_sequence(sequence, ref=lcs)
#         print(create_fasta(seq_id, shifted_sequence))
else:
    print("Either a reference sequence or a position must be specified.")

