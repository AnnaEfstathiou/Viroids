from Levenshtein import distance

def fasta_to_dict(fasta_file):
    '''
    Create a dictionary from a fasta file, where the headers are the keys and the values the sequences.
    '''
    sequences = {}
    current_id = None
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):  # Check if the line is a sequence identifier
            current_id = line[1:]  # Extract the identifier (excluding the '>' character)
            sequences[current_id] = []  # Initialize a new entry in the dictionary
        else:
            sequences[current_id].append(line)  # Append the sequence line to the current identifier's list
    # Join all parts of sequences together and return the dictionary
    return {seq_id: ''.join(seq_parts) for seq_id, seq_parts in sequences.items()}


def create_fasta(seq_id, sequence, line_width=60):
    '''
    Convert a sequence into FASTA format.
    '''
    return f">{seq_id}\n" + "\n".join(sequence[i:i+line_width] for i in range(0, len(sequence), line_width))


def reverse_complement(seq):
    '''
    Calculate the reverse complement of a DNA sequence.
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} # Mapping of complements
    return ''.join(complement.get(base, base) for base in reversed(seq)) # Return the reverse complement


def find_closest_match(seq, ref, max_dist=3):
    '''
    Find the closest match with fuzzy matching.
    '''
    best_pos = None
    best_dist = max_dist + 1  # Initialize best_dist to be just outside the acceptable range
    for i in range(len(seq) - len(ref) + 1):  # Iterate through each possible starting position in seq
        segment = seq[i:i+len(ref)]  # Extract the segment of seq that is the same length as ref
        dist = distance(segment, ref)  # Calculate the Levenshtein distance between the segment and ref
        if dist < best_dist:  # If this distance is the smallest so far
            best_dist = dist
            best_pos = i
    return best_pos if best_dist <= max_dist else -1  # Return the start position of the best match, or -1 if no suitable match was found

def shift_sequence(seq, pos=None, ref=None, max_dist = None):
    '''
    Shift the sequence either based on a given position or a reference sequence.
    '''
    if pos:  # If a position is specified
        return seq[pos-1:] + seq[:pos-1]  # Rotate the sequence to start from the specified position
    elif ref:  # If a reference sequence is specified
        pos = seq.find(ref)  # Find the position of the reference in the sequence
        # pos = find_closest_match(seq, ref, max_dist=3)
        if pos != -1:
            return seq[pos:] + seq[:pos]  # Rotate the sequence to align with the reference
        else:
            ref_rc = reverse_complement(ref)  # Get the reverse complement of the reference
            pos = seq.find(ref_rc)  # Find the position of the reverse complement in the sequence
            # pos = find_closest_match(seq, ref, max_dist=3)
            if pos != -1:
                return seq[pos:] + seq[:pos]  # Rotate the sequence to align with the reverse complement
            else:
                pos = find_closest_match(seq, ref, max_dist)
                return seq[pos:] + seq[:pos]
    else:
        raise ValueError("Either position or reference must be specified.")  # Raise an error if neither is specified
