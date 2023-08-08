import csv
import time
import numpy as np
import datetime as dt
from aligner import aligner
from constant_cord import constant_coord
def mismatch_filter(working_folder_name, to_filter, constant_region_loc,
                    allowed_mismatch_counts, expected_sequence, stats_file_address):
    # SETUP
    constant_region_numbers = list(np.arange(len(allowed_mismatch_counts), dtype=int))
    # Count sequences
    read_count = 0
    pass_count = 0

    # Holders for sequences to be written to txt files
    passed_sequence_holder = []
    failed_sequence_holder = []

    # Input file
    input_file_address = working_folder_name + '/' + to_filter

    # Output file names
    output_file_name = str(dt.date.today()) + ' ' + to_filter + ' mismatch filtered.txt'
    output_file_address = working_folder_name + '//' + output_file_name
    failed_output_file_name = str(dt.date.today()) + ' ' + to_filter + ' mismatch failed.txt'
    failed_output_file_address = working_folder_name + '//' + failed_output_file_name


    print('Started at', dt.datetime.now().strftime('%I:%M:%S'))
    start_time = time.time()


    # Open input and output files
    with open(output_file_address, 'a', newline='\n') as output_file, \
    open(failed_output_file_address, 'a', newline='\n') as failed_output_file:

        # Clear output files from any previous runs
        output_file.truncate(0)
        failed_output_file.truncate(0)

    # Align each line in input file with the expected sequence
        alignment_objects=aligner(input_file_address,expected_sequence)

    # Iterate through alignment objects and check for mismatches
        for alignment in alignment_objects:
            read_count+=1
            for locs, max_mismatches, i in zip(constant_region_loc, allowed_mismatch_counts,constant_region_numbers):
                mismatch_pass = True
                mismatches = 0
                cons_coord = []
                for loc in locs:
                    cons_coord.append(constant_coord(alignment[0],loc))
                partial_expected_sequence=alignment[0][cons_coord[0]:(cons_coord[1]+1)]
                partial_query_sequence=alignment[1][cons_coord[0]:(cons_coord[1]+1)]
                if partial_expected_sequence != partial_query_sequence:
                    for c, e in zip(partial_query_sequence,partial_expected_sequence):
                        if c != e:
                            mismatches += 1
                            if mismatches > max_mismatches:
                                mismatch_pass = False
                                failed_sequence_holder.append(alignment.query)
                                break
            if mismatch_pass:
                pass_count += 1
                # Add to holder
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

        # Write all sequences left in holders at the end of the loop
        # Passing sequence holder
        for seq in passed_sequence_holder:
            output_file.write(seq)
        passed_sequence_holder.clear()
        print('{:,} sequences passed, {:,} failed'.format(pass_count,
              read_count-pass_count))
        for seq in failed_sequence_holder:
            failed_output_file.write(seq)
        failed_sequence_holder.clear()
        
        end_time = time.time()
    print('\n{:,} reads processed in {:.0f} seconds'.format(read_count, end_time-start_time))
    print('{:,} sequences passed ({:.2%})'.format(pass_count, pass_count/read_count))
    print('Written to ', output_file_name)

    # Write the results to a .csv file
    # Opening the file in 'a' (append) mode means that lines are added to the end of the
    # existing file.

    # The structure of this file (set up in Main) is:
    # ['Biosample', 'Pct pass', 'Total reads', 'Passing reads', 'Region 1 failed'...]

    # Compile filtering statistics into stats
    stats = [to_filter, (pass_count/read_count), read_count, pass_count]
    # Write to file as a new line
    with open(stats_file_address, 'a', newline='') as stats_file:
        stats_writer = csv.writer(stats_file)
        stats_writer.writerow(stats)

    return output_file_name
