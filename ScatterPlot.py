import matplotlib.pyplot as plt
import numpy as np
import csv
def to_data(x_file, y_file):
    with open(x_file) as X, open(y_file) as Y:
        x_reader = csv.reader(X)
        y_reader = csv.reader(Y)
        next(x_reader)
        next(y_reader)
        x_count = 0
        y_count = 0
        for row in x_reader:
            x_count += int(row[1])
        for row in y_reader: 
            y_count += int(row[1])
        counts = [x_count, y_count]
        print(counts)
    with open(x_file) as X, open(y_file) as Y:
        x_reader = csv.reader(X)
        y_reader = csv.reader(Y)
        readers = [x_reader,y_reader]
        x_holder = {}
        y_holder = {}
        holders = [x_holder, y_holder]
        x_data = []
        y_data = []
        next(x_reader)
        next(y_reader)
        for reader, holder in zip(readers, holders):
            for row in reader:
                a=int(row[1])
                holder[row[0]]=a
        for x_key in x_holder:
            for y_key in y_holder:
                if x_key == y_key:
                    x_data.append(x_holder[x_key]/(x_count/1e6))
                    y_data.append(y_holder[y_key]/(y_count/1e6))
    return(x_data,y_data)

def scatter_plot(x_data, y_data,x_label,y_label):

    # Convert data to numpy arrays
    x_data = np.array(x_data)
    y_data = np.array(y_data)

    fig, ax = plt.subplots()

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.scatter(x_data,y_data)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    plt.show

    # Show the plot
    plt.show()
x_data_path='/Users/veritas/Desktop/HYR/Analysis/L2/Results/L2S5_R1.fastq Length filtered.txtCount possible sequence stats.csv'
y_data_path='/Users/veritas/Desktop/HYR/Analysis/L2/Results/L2S6_R1.fastq Length filtered.txtCount possible sequence stats.csv'
data = to_data(x_data_path,y_data_path)
scatter_plot(data[0],data[1],'L2S5','L2S6')





