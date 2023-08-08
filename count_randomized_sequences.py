import numpy as np
import csv
from Bio import Align,AlignIO
from constant_cord import constant_coord

def count_randomized_sequences(library, to_counter, working_folder_name,
                               expected_sequence):
    """Count the frequency of each randomized sequence in the filtered data for each
    biosample.

    Parameters
    ---
    biosamples : list
        Biosample names
    filtered_file_list : list
        Filtered files to read
    working_folder_name : str
        Working folder
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [0, 0...]]
    randomized_bases : list
        Randomized base coordinates

    Returns
    ---
    possible_sequence_dict : dict
        Dictionary of the form ['Randomized sequence': [Lib_count, sel_1_count...]]
    read_count : list
        List of the total identified reads for each biosample.
    stats : list
        Stats summarizing the filtering results
        [[Biosample, processed_count, found, not_found]...]

    Notes
    ---
    Reads filtered files and updates possible_sequence_dict with counts

    Using randomized sequences as dictionary keys is much faster than using a list
    because dictionary keys are hashed, making searching much more efficient.
    """

    # Setup
    # stats holds statistics describing how many read matched expected library sequences
    # for each biosample. [[biosample, processed_count, found, not_found]...]
    stats = []
    # read_counts holds the number of identified reads which matched expected library
    # sequences for each biosample. This will be used later to calculate enrichment
    # factors.
    read_counts = []

    # Iterate through biosamples
        # Track how many sequences were found in the expected library sequences

    aligner=Align.PairwiseAligner()
    aligner.mismatch_score=0
    aligner.open_gap_score=-1
    aligner.extend_gap_score=0
    aligner.match_score=1


    print('Processing: ' + to_counter)
    random_seq_holder = {}
    to_file=[]
    processed_count = 0
    # Open the filtered sequence file and read each sequence.
    filtered_file_address = working_folder_name + '/' + to_counter
    with open(filtered_file_address) as filtered_file:
        for sequence in filtered_file:
            processed_count += 1
            # Extract randomized bases
            extract_random = ''
            if library == 'L5':
                alignments = aligner.align(expected_sequence,sequence)
                alignment=alignments[0]
                constant_column1 = constant_coord(alignment[0],123)
                constant_column2 = constant_coord(alignment[0],133)
                a=alignment[1][:constant_column1+1]
                b=alignment[1][:constant_column2+7]
                constant_region=alignment[1][constant_column1:(constant_column2+1)]
                query_sequence1 = a.replace('-','')
                query_sequence2 = b.replace('-','')
                extract_random=(query_sequence1[-6:-1]+ constant_region + query_sequence2[-6:-1])
                extract_random = extract_random.replace('-','')
            # Exclude sequences which are not in the expected sequence list from analysis
            elif library == 'L2':
                alignments = aligner.align(expected_sequence,sequence)
                alignment=alignments[0]
                constant_column1 = constant_coord(alignment[0],117)
                constant_column2 = constant_coord(alignment[0],139)
                extract_random = alignment[1][constant_column1:constant_column2+1]
                extract_random = extract_random.replace('-','')
            if extract_random in random_seq_holder.keys():
                random_seq_holder[extract_random]+=1
            else: 
                random_seq_holder[extract_random]=1
    # Counting statistics
    for key in random_seq_holder:
        to_file.append([key,random_seq_holder[key]])
    print('{:,} sequences processed'.format(processed_count))
    # Write counting statistics to a .csv file
    headers = ['sequence','occurence']
    output_file_name = to_counter + 'Count possible sequence stats.csv'
    output_file_address = working_folder_name + '/Results/' + output_file_name
    with open(output_file_address, 'w', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(headers)
        for count in to_file:
            stats_writer.writerow(count)



