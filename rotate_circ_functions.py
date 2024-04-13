from Levenshtein import distance

def fasta_to_dict(fasta_file):
    '''
    Create a dictionary from a fasta file, where the headers are the keys and the values the sequences.
    '''
    sequences = {}
    current_id = None
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):  # check if the line is a sequence header
            current_id = line[1:]  # extract the header (excluding the '>' character)
            sequences[current_id] = []  # initialize a new entry in the dictionary
        else:
            sequences[current_id].append(line)  # append the sequence line to the current identifier's list
    # join all parts of sequences together and return the dictionary
    return {seq_id: ''.join(seq_parts) for seq_id, seq_parts in sequences.items()}


def reverse_complement(seq):
    '''
    Calculate the reverse complement of a DNA sequence.
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} # mapping of complements
    return ''.join(complement.get(base, base) for base in reversed(seq)) # return the reverse complement


def find_closest_match(seq, ref, max_dist=3):
    '''
    Find the closest match with fuzzy matching.
    Levenshtein distance: string metric for measuring the difference between two sequences. 
    The distance between them is the minimum number of single-character edits (insertions, deletions or substitutions) required to change one sequence into the other.
    '''
    best_pos = None
    best_dist = max_dist + 1  # initialize best_dist to be just outside the acceptable range
    for i in range(len(seq) - len(ref) + 1):  # iterate through each possible starting position in sequence
        segment = seq[i:i+len(ref)]  # extract the segment of seq that is the same length as reference
        dist = distance(segment, ref)  # calculate the Levenshtein distance between the segment and reference
        if dist < best_dist:  # if this distance is the smallest so far
            best_dist = dist
            best_pos = i
    return best_pos if best_dist <= max_dist else -1  # return the start position of the best match, or -1 if no suitable match was found


def shift_sequence(seq, pos=None, ref=None, max_dist = None):
    '''
    Shift the sequence either based on a given position or a reference sequence.
    '''
    if pos:  # ff a position is specified
        return seq[pos-1:] + seq[:pos-1]  # rotate the sequence to start from the specified position
    elif ref:  # if a reference sequence is specified
        pos = seq.find(ref)  # find the position of the reference in the sequence
        if pos != -1: # if the position is found
            return seq[pos:] + seq[:pos]  # rotate the sequence to align with the reference
        else:
            ref_rc = reverse_complement(ref)  # get the reverse complement of the reference
            pos = seq.find(ref_rc)  # find the position of the reverse complement in the sequence
            if pos != -1: # if the position is found
                return seq[pos:] + seq[:pos]  # rotate the sequence to align with the reverse complement
            else: # if the position is not found
                pos = find_closest_match(seq, ref, max_dist) # find the closest position that coresponds to reference
                return seq[pos:] + seq[:pos] # rotate the sequence to align with the fuzzing matching
    else:
        raise ValueError("Either position or reference must be specified.")  


def write_fasta(fasta_file, seq_id, sequence):
    '''
    Write fasta format to file.
    '''
    fasta_file.write(f">{seq_id}\n{sequence}\n")