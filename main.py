import os
import datetime as dt
import csv
from mismatch_filter import mismatch_filter
from read_parameter import read_parameter_file
from identify_fastq import identify_fastq_files
from Q_score_filter import Q_score_filter
from Alignment_score_filter import Alignment_score_filter
from count_randomized_sequences import count_randomized_sequences
from Length_filter import Length_filter


# SETUP

# Enter parameter file name and path here
parameter_file_path = "/Users/veritas/Desktop/HYR/Analysis/aptamer_parameter_L2.csv"

# Funtion on/off switches

# find_fastq_files identifies the files associated with each biosample
# Always set to True
run_find_fastq_files = True
run_length_filter = True
# Filtering
# Q_score_filter and mismatch_filter read an input file and filter by data quality
# If run_Q_score_filter is set to False but run_mismatch_filter is set to True, the
# program looks for previously Q score filtered text files and uses those as the input
# for mismatch filtering
# If both are set to False but a subsequent step is set to True, the program looks for
# previously mismatch filtered text files and uses those as the input for future steps
run_Q_score_filter = False
run_alignment_score_filter = False
run_mismatch_filter = False

#
run_count_sequences = False




# 1. READ IN PARAMETERS
# read_parameter_file returns p : {'parameter_name': parameter_value}
print('\n---\nReading parameters\n---')
p=read_parameter_file(parameter_file_path)

# Verify and/or create required folders
required_folders = ['Plots', 'Results'] #result folder, plot folder
for folder in required_folders:
    full_path = p['working_folder_name'] + '//' + folder
    if os.path.exists(full_path) == False:
        os.makedirs(full_path)

# 2. IDENTIFY ALL .FASTQ FILES
# # Returns fastq_files: [[biosample_1_file_1, biosample_1_file_2...],
# [biosample_2_file_1, biosample_2_file_2...]...]
# This accomodates cases where one biosample is associated with multiple .fastq files.
if run_find_fastq_files:
    print('\n---\nLooking for .fastq files in:\n' + p['working_folder_name'])
    fastq_files = identify_fastq_files(p['working_folder_name'], p['lib_name'])

# Length filter

if run_length_filter:
    length_filtered_files = []
    stats_file_name = str(dt.date.today()) + ' Length filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '/Results/' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        header_writer = csv.writer(stats_file)
        headers = ['Library', 'Round' 'Pct pass', 'Total reads', 'Passing reads',
                    'Fail count']
        header_writer.writerow(headers)
    for to_filter in fastq_files:
        print('\nLength filtering:\n' + p['lib_name'])
        filtered_file_name = Length_filter(p['working_folder_name'], p['lib_name'], to_filter,
                                            stats_file_address)
        length_filtered_files.append(filtered_file_name)
    for to_filter in length_filtered_files:
        # Count in each mismatch-filtered sequence file
        print('\n---\nCounting randomized sequences\n---')
        count_randomized_sequences(p['lib_name'],to_filter,p['working_folder_name']
                                , p['expected_sequence'])

# 3. QUALITY FILTER
if run_Q_score_filter:
    print('\n---\nQuality filtering\n---')
    print('\nQ score filtering parameters for the entire sequence')
    print('\nNo more than {} bases below Q{}'.format(p['QF'], p['Q1']))
    print('and no bases below Q{}'.format(p['Q0']))

    # Set up a list of quality-filtered files
    Q_score_filtered_files = []

    # Set up a .csv file to hold the Q score filtering stats
    stats_file_name = str(dt.date.today()) + ' Q score filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '/Results/' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        header_writer = csv.writer(stats_file)
        headers = ['Library', 'Round' 'Pct pass', 'Total reads', 'Passing reads',
                    'Q0 fail', 'Q1 fail']
        header_writer.writerow(headers)

    # Call Q_score_filter for each biosample and its associated .fastq file(s)
    # Returns the name of the Q score filtered file and writes stats to the .csv file
    for to_filter in fastq_files:
        print('\nQ score filtering:\n' + p['lib_name'])
        filtered_file_name = Q_score_filter(p['working_folder_name'], p['lib_name'], to_filter,
                                            p['Q1'], p['QF'], p['Q0'],
                                            stats_file_address)
        Q_score_filtered_files.append(filtered_file_name)

# 4. ALIGNMENT SCORE FILTER

if run_alignment_score_filter:
    print('\n---\nAlignment filtering\n---')
    print('\nAlignment score filtering parameters')
    print('\nNo sequence with alignment score below {}'.format(p['A0']))

    # Set up a list of alignment-filtered files
    Alignment_filtered_files = []

    # Set up a .csv file to hold the Alignment filtering stats
    stats_file_name = str(dt.date.today()) + ' Alignment filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '/Results/' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        header_writer = csv.writer(stats_file)
        headers = ['Library', 'round', 'Pct pass', 'Total reads', 'Passing reads']
        header_writer.writerow(headers)

    # Call Alignment_score_filter for each .fastq file(s)
    # Returns the name of the Alignmnet score filtered file and writes stats to the .csv file
    for to_filter in Q_score_filtered_files:
        print('\nAlignment score filtering:\n' + p['lib_name'])
        filtered_file_name = Alignment_score_filter(p['working_folder_name'],p['lib_name'],to_filter,p['A0'],p['expected_sequence'],stats_file_address)
        Alignment_filtered_files.append(filtered_file_name)

# 5. MISMATCH FILTER
if run_mismatch_filter:
    print('\n---\nMismatch filtering\n---')
    print('Filters set to:')
    for m in range(len(p['allowed_mismatch_counts'])):
        print('> Constant region {}: allow {} mismatches'.format(
                m+1, p['allowed_mismatch_counts'][m]))

    # Set up a list of mismatch-filtered files
    mismatch_filtered_files = []

    # Set up a .csv file to hold the Q score filtering stats
    stats_file_name = str(dt.date.today()) + ' mismatch filtering stats.csv'
    stats_file_address = p['working_folder_name'] + '/Results/' + stats_file_name
    with open(stats_file_address, 'w', newline='') as stats_file:
        # Clear the whole file
        stats_file.truncate(0)
        # Write the headers
        header_writer = csv.writer(stats_file)
        headers = ['Round', 'Pct pass', 'Total reads', 'Passing reads']
        header_writer.writerow(headers)

    # Call mismatch_filter for each biosample and its associated .fastq file(s)
    # Returns the name of the mismatch filtered file and stats describing how many
    # sequences passed.
    # mismatch_filter also automatically generates plots showing base and mismatch
    # frequency at each position and writes the same data to .csv files.
    for to_filter in Alignment_filtered_files:
        print('\nMismatch filtering:\n'+ to_filter)
        filtered_file_name = mismatch_filter(p['working_folder_name'],to_filter,
                                            p['constant_region_loc'],p['allowed_mismatch_counts'],p['expected_sequence'],
                                            stats_file_address)
        mismatch_filtered_files.append(filtered_file_name)

#6 COUNT RANDOMIZED SEQUENCE
if run_count_sequences:
    for to_filter in mismatch_filtered_files:
        # Count in each mismatch-filtered sequence file
        print('\n---\nCounting randomized sequences\n---')
        count_randomized_sequences(p['lib_name'],to_filter,p['working_folder_name']
                                , p['expected_sequence'])

