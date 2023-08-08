import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from identify_csv_files import identify_csv_files
def toFasta(working_folder_path, library_name, top_number, repetition):
    files=identify_csv_files(working_folder_path,library_name)
    for file in files:
        trial_name_index = file.find(library_name + 'S')
        trial_name = file[trial_name_index : trial_name_index+5]
        with open(working_folder_path +'/' + file,'r') as csv_file:
            dict={}
            reader=csv.reader(csv_file)
            next(reader)
            for row in reader:
                dict[row[0]]=int(row[1])
            sorted_dict = sorted(dict.items(), key=lambda x:x[1], reverse=True)
            if top_number != 0:
                tops = sorted_dict[0:top_number]
            else: 
                tops = sorted_dict
            if repetition:
                record_list = []
                index = 1
                for item in tops:
                    record = SeqRecord(Seq(str(item[0])),id = str(index),description = trial_name + 'top' + str(index))
                    index += 1
                    for i in range(0,int(item[1])):
                        record_list.append(record)
                        i+=1
                SeqIO.write(record_list,working_folder_path + '/' + trial_name + 'top' + str(top_number) + '.fasta', 'fasta')
            else:
                record_list = [] 
                index = 1
                for item in tops:
                    record = SeqRecord(Seq(str(item[0])),id = str(index), description = trial_name + 'top' + str(index))
                    index += 1
                    record_list.append(record)
                SeqIO.write(record_list,working_folder_path + '/' + trial_name + 'top' + str(top_number) + '.fasta', 'fasta')

toFasta('/Users/veritas/Desktop/HYR/Analysis/L2/Results', 'L2', 0, True)
