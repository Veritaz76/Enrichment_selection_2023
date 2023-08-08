import csv
import datetime as dt
def read_parameter_file(parameter_file_path):
    """Read in data processing parameters from a .csv file and output to dictionary

    Parameters
    ---
    parameter_file_name : str
        full path and file name of parameters.csv

    Returns
    ---
    p : a dictionary of parameter names and values
        {parameter_name: parameter_value}

    Notes
    ---
    Converts base index values from SeqBuilder-type counting, which starts at 1,
    to Python-type counting, which starts at 0. (In parameter .csv file, enter sequence counting starting at 1)

    The input parameters.csv file is expected to have the following format
    parameter_name, description, value1, value2...

    The parameters used will be written to output .csv files as well for record-keeping.

    Blank lines are okay, but any line with an entry in the parameter_name (first)
    column is assumed to contain a parameter to be read.
    """
    print('Reading parameters from '+ parameter_file_path)

    # Setting up an empty dictionary to hold lines read from .csv file
    # After reading the parameter file p0 will look like {parameter_name: [value1, value2...]}
    p0 = {}

    # 1. Read parameter file
    # - Tracking all rows and write to another .csv file for record-keeping
    rows=[]
    with open(parameter_file_path,'r') as parameter_file:
        parameter_file_reader=csv.reader(parameter_file) #Initialize scanner
        for row in parameter_file_reader:
            rows.append(row)
            # Only read rows with a value in the parameter_name colunm
            if row[0] != '':
                values = []
                for value in row[2:]:
                    if value != '':
                        values.append(value)
                if len (values) == 0:
                    # Blanked values are default to be None
                    p0[row[0]] = None
                else: 
                    p0[row[0]] = values
    p={} # Initializing parameter dictionary for output

    # Parameter type 1: simple parameters (sp)
    # Blank values default to None
    for sp in ['lib_name','working_folder_name', 'expected_sequence',
               'twist_sequence_file', 'output_format', 'full_seq_format']:
                if (sp in p0.keys() and p0[sp]):
                    p[sp] = p0[sp][0]
                else: 
                    p[sp] = None
    # Parameter type 2: simple numerical parameters
    # Blank values default to 0
    for np in ['Q0', 'Q1', 'Q2', 'QF','A0','exclude_at_start', 'exclude_at_end', 'min_lib_count']:
        if (np in p0.keys() and p0[np]):
            p[np] = int(p0[np][0])
        else: 
            p[np] = 0
    # Correct Q2 to default to Q1 so failed bases are not allowed in randomized regions        
    if p['Q2'] == 0:
        p['Q2'] == p['Q1']
    # Correct min_lib_count to 1 to avoid divide by zero errors later
    if p['min_lib_count'] == 0:
        p['min_lib_count'] == 0    
    

    # Numerical lists
    p['allowed_mismatch_counts'] = []
    for m in p0['allowed_mismatch_counts']:
        p['allowed_mismatch_counts'].append(int(m))
    # Coordinate lists which need integer extraction and counting adjustment
    # tRNA and constant region coordinates (tRNA coords are optional)
    if 'tRNA_coords' in p0:
        p['tRNA_coords'] = [int(p0['tRNA_start'][0])-1, int(p0['tRNA_end'][0])]
    else:
        p['tRNA_coords'] = None
    p['constant_region_loc'] = []
    for i, j in zip(p0['constant_region_starts'], p0['constant_region_ends']):
        p['constant_region_loc'].append([int(i)-1, int(j)])
    # List of all randomized bases
    p['randomized_bases'] = []
    for b in p0['randomized_bases']:
        p['randomized_bases'].append(int(b)-1)
    # Expected pairs in terms of sequence coordinates
    if (p0['pairs_5prime'] and p0['pairs_3prime']):
        p['expected_pairs'] = []
        for a, b in zip(p0['pairs_5prime'], p0['pairs_3prime']):
            p['expected_pairs'].append([int(a)-1,int(b)-1])
    else:
        p['expected_pairs'] = None
    # Expected pairs in terms of randomized base coordinates
    p['rand_expected_pairs'] = []
    if p['expected_pairs']:
        for c in p['expected_pairs']:
            p['rand_expected_pairs'].append([p['randomized_bases'].index(c[0]),
                                            p['randomized_bases'].index(c[1])])

    # Print output
    print('\nParameters set to:')
    for parameter in p:
        print('> {} : {}'.format(parameter, p[parameter]))

    # Write the initial and adjusted parameters used to .csv files in the working folder
    # for record-keeping purposes
    # date lib_name parameters as entered.csv will be an exact copy of the input file
    as_entered_file_name = str(dt.date.today()) + ' ' + p['lib_name'] + ' parameters as entered.csv'
    as_entered_file_address = p['working_folder_name'] + '//' + as_entered_file_name
    with open(as_entered_file_address, 'w', newline='') as output_file:
        parameter_writer = csv.writer(output_file)
        for row in rows:
            parameter_writer.writerow(row)

    # date lib_name parameters adjusted.csv will be the parameters after adjstment
    adjusted_file_name = str(dt.date.today()) + ' ' + p['lib_name'] + ' parameters adjusted.csv'
    adjusted_file_address = p['working_folder_name'] + '//' + adjusted_file_name
    with open(adjusted_file_address, 'w', newline='') as output_file:
        parameter_writer = csv.writer(output_file)
        for parameter_name in p:
            parameter_writer.writerow([parameter_name, p[parameter_name]])

    return p