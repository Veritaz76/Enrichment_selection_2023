import numpy as np
from Bio import Align,AlignIO
def aligner(input_file_address, expected_sequence):

    # Holders for sequences to be written to txt files
    alignment_objects_holder = []

    # Setup Aligner
    aligner=Align.PairwiseAligner()
    target=expected_sequence
    aligner.mismatch_score=0
    aligner.open_gap_score=-2
    aligner.extend_gap_score=-1
    aligner.match_score=1

    # Open input file
    # Pairwise align and append to alignment_object_holder
    with open(input_file_address, 'r', newline='\n') as input_file:
        for query in input_file:
            alignments=aligner.align(target,query)
            alignment=alignments[0]
            alignment_objects_holder.append(alignment)
    
    # Return all alignment objects

    return alignment_objects_holder

