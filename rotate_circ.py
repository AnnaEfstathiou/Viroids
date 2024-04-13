import argparse
import os

def main(args):

    ''' CHECKING FOR ERRORS '''

    # check for the presence of rotate_circ_functions.py in the current directory
    if not os.path.exists('rotate_circ_functions.py'):
        raise ImportError("rotate_circ_functions.py is not found in the current directory. Please ensure it is available.")

    # import the functions from rotate_circ_functions.py
    from rotate_circ_functions import fasta_to_dict, shift_sequence, write_fasta

    # check for the presence of the input file in the current directory
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"The input file {args.input} does not exist.")
    
    if args.ref and not os.path.exists(args.ref):
        raise FileNotFoundError(f"The reference file {args.ref} does not exist.")

    ''' MAIN CODE '''

    with open(args.input, 'r') as input_file:
        fasta_sequences = fasta_to_dict(input_file)

    with open(args.output, 'w') as output_file:
        # Process sequences based on the reference or position
        if args.ref:
            with open(args.ref, 'r') as ref_file:
                reference_sequences = fasta_to_dict(ref_file)
                reference_sequence = next(iter(reference_sequences.values()))
                for seq_id, sequence in fasta_sequences.items():
                    try:
                        shifted_sequence = shift_sequence(sequence, ref=reference_sequence, max_dist=args.max_dist)
                        write_fasta(output_file, seq_id, shifted_sequence)
                    except ValueError as e:
                        print(f"Error for {seq_id}: {str(e)}")
        elif args.pos:
            for seq_id, sequence in fasta_sequences.items():
                shifted_sequence = shift_sequence(sequence, pos=args.pos)
                write_fasta(output_file, seq_id, shifted_sequence)
        else:
            raise ValueError("Either a reference sequence or a position must be specified.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Shift circular FASTA sequences based on a reference or a position.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA file.")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output FASTA file.")
    parser.add_argument("-r", "--ref", type=str, help="Reference FASTA file.")
    parser.add_argument("-p", "--pos", type=int, help="Position to shift the sequence.")
    parser.add_argument("-d", "--max_dist", type=int, default=3, help="Maximum distance of the Levenshtein distance. The default value is 3.")
    args = parser.parse_args()
    main(args)