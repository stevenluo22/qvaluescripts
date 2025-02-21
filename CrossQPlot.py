#Data collection and analysis part
import prody as pr
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.colors import LinearSegmentedColormap

import csv

def run(args):
    # Define the file name
    args.inputCSV

    # Initialize the dictionary
    cross_q_val_table = {}

    # Read the CSV file
    with open(file_name, mode='r') as file:
        reader = csv.reader(file)
        headers = next(reader)  # Skip the header row

        # Populate the dictionary
        for row in reader:
            row_key = int(row[0])  # Use the first column as the key
            cross_q_val_table[row_key] = {}

            # Iterate over the remaining items in the row
            for col_index, value in enumerate(row[1:], start=1):
                cross_q_val_table[row_key][int(headers[col_index])] = float(value)

    import numpy as np

    # Convert dictionary to a 2D NumPy array
    #q_val_array = np.array([[cross_q_val_table[i][j] for j in range(1, 9)] for i in range(1, 9)])

    # Define custom colormap going from red to yellow to green
    colors = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]  # Red to Yellow to Green
    cmap_name = 'red_yellow_green'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

    file_name = args.outputplot
    q_val_map = sns.clustermap(cross_q_val_table, cmap = cm)

    # Access the colorbar and set its label
    cbar = q_val_map.ax_heatmap.collections[0].colorbar
    cbar.set_label('Mutual Q value', fontsize=15)

    # Axes
    ax = q_val_map.ax_heatmap
    #plt.title("Mutual Q values, Mixed Memory (weights single memories 100), 500 K to 200 K")
    # Access the data array used in the clustermap
    data = q_val_map.data2d.values

    # Loop through the data array and annotate each cell with its value
    for i in range(len(cross_q_val_table)):
        for j in range(len(cross_q_val_table)):
            ax.text(j + 0.5, i + 0.5, '{:.2f}'.format(data[i, j]), ha='center', va='center', color='black', fontsize=args.fontSize)
    plt.savefig(file_name)

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Draw plots from cross Q data csv")
    parser.add_argument("-i", "--inputCSV", help="csv file with list of pdb files paths", default="CrossQ.csv", type=str)
    parser.add_argument("-o", "--outputplot", help="Name of output plot", default="CrossQ.jpg", type=str)
    parser.add_argument("-f", "--fontSize", help="font size", default=10, type=float)

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    run(args)
    
if __name__=="__main__":
    main()