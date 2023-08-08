import os
import csv
import time
import datetime as dt
from aligner import aligner
from Bio import Align,AlignIO
def Alignment_score_filter(working_folder_name, library, to_filter, A0,
                         expected_sequence,stats_file_address):
# Setup

    # Count sequences which fall into each category
    pass_count = 0
    fail_count = 0
    read_count = 0

    # Holders for sequences to be written to .txt files
    passed_sequence_holder = []
    failed_sequence_holder = []
    
    start_time = time.time()

    # Input file
    input_file_address = working_folder_name + '/' + to_filter

    # File names and addresses
    output_file_name = str(dt.date.today()) + ' ' + library + ' ' + to_filter[13:19] + ' Alignment score filtered.txt'
    output_file_address = working_folder_name + '/' + output_file_name
    failed_output_file_name = str(dt.date.today()) + ' ' + library + to_filter + ' Alignment score failed.txt'
    failed_output_file_address = working_folder_name + '/' + failed_output_file_name

    print('Started at', dt.datetime.now().strftime('%I:%M:%S'))
    start_time = time.time()

    # Open input and output files
    with open(output_file_address, 'a', newline='\n') as output_file, \
         open(failed_output_file_address, 'a', newline='\n') as failed_output_file:

        # Clear output files from any previous runs
        output_file.truncate(0)
        failed_output_file.truncate(0)

        print ('\nProcessing file:\n' + to_filter)
        print ('Started filtering reads at',dt.datetime.now().strftime('%I:%M:%S'))

        # Align
        alignment_objects=aligner(input_file_address,expected_sequence)
        

        # Filter 
        for alignment in alignment_objects:
            read_count+=1
            if alignment.score < A0:
                fail_count+=1
                failed_sequence_holder.append(alignment.query)
            else:
                pass_count+=1
                passed_sequence_holder.append(alignment.query)
            
            # Whenever any holder is full, write the sequences in the holder to the
            # appropriate file and clear the holder.
            # This is done to avoid generating enormous lists which exceed the
            # computer's available memory.
            holder_length = 100000

            if len(passed_sequence_holder) >= holder_length:
                for seq in passed_sequence_holder:
                    output_file.write(seq)
                passed_sequence_holder.clear()
                print('{:,} sequences passed, {:,} failed'.format(pass_count,
                      read_count-pass_count))

            elif len(failed_sequence_holder) >= holder_length:
                for seq in failed_sequence_holder:
                    failed_output_file.write(seq)
                failed_sequence_holder.clear()
        
        # At the end of the loop, write the sequences in each holder to the
        # appropriate file
        for seq in passed_sequence_holder:
            output_file.write(seq)
        passed_sequence_holder.clear()
        # print('{:,} sequences passed, {:,} failed ({:.2% pass})'.format(pass_count,
        #      read_count-pass_count, pass_count/read_count))

        for seq in failed_sequence_holder:
            failed_output_file.write(seq)
        failed_sequence_holder.clear()

    end_time = time.time()
    
    print ('\nAlignment score filtering results for ' + library + to_filter)
    print ('{:,} reads processed'.format(read_count))
    print('{:,} sequences failed A0 = {} ({:.2%})'.format(fail_count, A0,
          fail_count/read_count))
    print('\nFiltered sequences written to:\n' + output_file_name)
    print('\nTotal run time = {:.0f} seconds ({:.1f} minutes)'.format(end_time-start_time,
          (end_time-start_time)/60))        
            
    # Write the results to a .csv file
    # Opening the file in 'a' (append) mode means that lines are added to the end of the
    # existing file.

    # The structure of this file (set up in Main) is:
    # ['Biosample', 'Pct pass', 'Total reads', 'Passing reads',
    # 'Short sequences', 'Q0 fail', 'Q1 fail', 'Q2 fail']

    # Compile filtering statistics into stats
    stats = [library, to_filter, (pass_count/read_count), read_count, pass_count]
    with open(stats_file_address, 'a', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(stats)

    return output_file_name