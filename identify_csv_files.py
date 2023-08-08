import os
def identify_csv_files(working_folder_path,library_name):
    found_file_names = os.listdir(working_folder_path)
    csv_files = []
    print('\nIdentified the following files:')
    end_name = '.csv'

    for f in found_file_names:
        if (library_name in f and f[-len(end_name):] == end_name and 'Count possible sequence' in f):
            print('>',f)
            csv_files.append(f)
    return sorted(csv_files,reverse = False)