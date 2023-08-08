import os

def identify_fastq_files(working_folder_path,library_name):
    """Finds all fastq files associated with each biosample in a list

    Parameters
    ---
    working_folder_path : str
        full path of the folder in which the files will be found
    biosamples : list
        list of biosamples for which to identify files

    Returns
    ---
    fastq_files : list
        [[biosample_1_file_1, biosample_1_file_2...],[biosample_2_file_1, biosample_2_file_2...]]

    Notes
    ---
    Set up to allow for one biosample having multiple associated .fastq files as is
    true for some NextSeq samples
    """
    # identify all files in working_folder_name and initiate a list of .fastq files
    found_file_names = os.listdir(working_folder_path)
    fastq_files = []
    print('\nIdentified the following files:')
    end_name = '.fastq'
    for f in found_file_names:
        if (library_name in f and f[-len(end_name):] == end_name):
            print('>',f)
            fastq_files.append(f)

    return fastq_files
