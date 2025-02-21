#Data collection and analysis part
import prody as pr
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import argparse
import importlib.util

def Distance(Vec_1, Vec_2):
    x1, y1, z1 = Vec_1
    x2, y2, z2 = Vec_2
    
    xcont = (x2 - x1)**2
    ycont = (y2 - y1)**2
    zcont = (z2 - z1)**2
    distance = np.sqrt(xcont + ycont + zcont)
    return distance

def TermQ(r_i, r_j, residue_i, residue_j, printTermQ=False):
    numerator = (r_i - r_j)**2/100
    Sigma_ij = 0.1*(abs(residue_i-residue_j)**0.15)
    if Sigma_ij == 0:
        return 0
    termQ = np.exp(-numerator/(2*Sigma_ij**2))
    if printTermQ:
        print(f"Between residues {residue_i} and {residue_j}, the term Q is {termQ}.")
    return termQ

def CrossQ(coords_a, coords_b, name_ids, chain_ids, residue_ids, printCrossQ=False, printTermQ=False):
    Contact_Cutoff = 3
    count_eligible_contact = 0
    TotalQ = 0
    for atom_i in range(len(coords_a)):
        for atom_j in range(len(coords_b)):
            if(abs(residue_ids[atom_i] - residue_ids[atom_j]) >= Contact_Cutoff and chain_ids[atom_i] == chain_ids[atom_j]):
                count_eligible_contact += 1
                Distance_a = Distance(coords_a[atom_i], coords_a[atom_j])
                Distance_b = Distance(coords_b[atom_i], coords_b[atom_j])
                TotalQ = TotalQ + TermQ(Distance_a, Distance_b, residue_ids[atom_i], residue_ids[atom_j], printTermQ)
    if printCrossQ:
        print(f"Amongst {count_eligible_contact}, we have a crossQ value of {TotalQ/count_eligible_contact}")
    return TotalQ/count_eligible_contact

def calculateQ(pathtoPDB1, pathtoPDB2, printCrossQ=False, printTermQ=False):
    structure_1 = pr.parsePDB(pathtoPDB1)

    # Extract coordinates, residue names, and atom names for the first structure
    coords_1 = structure_1.getCoords()  # numpy array of coordinates
    residues_1 = structure_1.getResnums()  # residue names
    atom_names_1 = structure_1.getNames()  # atom names
    chains_1 = structure_1.getChids()
    
    # Repeat for the second PDB file
    structure_2 = pr.parsePDB(pathtoPDB2)

    coords_2 = structure_2.getCoords()
    residues_2 = structure_2.getResnums()
    atom_names_2 = structure_2.getNames()
    chains_2 = structure_2.getChids()

    # Define a set of DNA residue names    
    # Initialize new lists for the filtered data
    CA_coords_1 = []
    CA_residues_1 = []
    CA_atom_names_1 = []
    CA_chains_1 = []
    CA_indicies_1 = []

    CA_coords_2 = []
    CA_residues_2 = []
    CA_atom_names_2 = []
    CA_chains_2 = []
    CA_indicies_2 = []

    #1 based indexing to align with fragment memory notations
    CA_index_1 = 1
    CA_index_2 = 1
    
    # Loop through the indices and filter based on the condition
    for i in range(len(atom_names_1)):
        if atom_names_1[i] == "CA":
            CA_coords_1.append(coords_1[i])
            CA_residues_1.append(residues_1[i])
            CA_atom_names_1.append(atom_names_1[i])
            CA_chains_1.append(chains_1[i])
            CA_indicies_1.append(CA_index_1)
            CA_index_1 += 1
    for i in range(len(atom_names_2)):
        if atom_names_2[i] == "CA":
            CA_coords_2.append(coords_2[i])
            CA_residues_2.append(residues_2[i])
            CA_atom_names_2.append(atom_names_2[i])
            CA_chains_2.append(chains_2[i])
            CA_indicies_2.append(CA_index_2)
            CA_index_2 += 1

    if (CA_index_1 != CA_index_2):
        print("Warning! Both files do not have the same number of atoms for Q values! Q Values may not be accurate")

    if printCrossQ:
        print(f"Calculating Cross Q values between {pathtoPDB1} and {pathtoPDB2}")
    return CrossQ(CA_coords_1, CA_coords_2, CA_atom_names_1, CA_chains_1, CA_residues_1, printCrossQ, printTermQ)

def run(args):
    spec = importlib.util.spec_from_file_location("fileList", args.pythonlist)
    fileList = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fileList)
    filelist = fileList.fileList()
    cross_q_val_table = {}
    for a in range(0, len(filelist)):
        cross_q_val_table[a] = {}
        for b in range(0, len(filelist)):
            cross_q_val_table[a][b] = calculateQ(filelist[a], filelist[b], args.printCrossQ, args.printTermQ)
            #print(filelist[a], filelist[b])
            
    # Define the file name
    file_name = args.outputCSV

    print(file_name)
    print(1/0)

    # Create a CSV writer object
    with open(file_name, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write the header row
        writer.writerow([''] + list(range(0, len(filelist))))  # Empty cell for the top-left corner

        # Write the data
        for i in range(0, len(filelist)):
            row = [i] + [cross_q_val_table[i][j] for j in range(0, len(filelist))]
            writer.writerow(row)

    print("cross_q_val_table saved to:", file_name)

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
        description="Calculating Cross-Q/Mutual-Q of pdb files")
    parser.add_argument("-y", "--pythonlist", help="List in form of python script", type=str)
    parser.add_argument("-p", "--pdblist", help="The name of the protein", type=str) #Not supported yet
    parser.add_argument("-t", "--txtlist", help="text file with list of pdb files paths", type=str) #Not supported yet
    parser.add_argument("-i", "--inputCSV", help="csv file with list of pdb files paths", type=str) #Not supported yet
    parser.add_argument("-c", "--outputCSV", help="Name of csv output data file", default="CrossQ.csv", type=str)
    parser.add_argument("-o", "--outputplot", help="Name of output plot", default="CrossQ.jpg", type=str)
    parser.add_argument("-f", "--fontSize", help="font size", default=13, type=float)
    parser.add_argument("--printTermQ", action="store_true", help="additional print each term Q", default=False)
    parser.add_argument("--printCrossQ", action="store_true", help="additional print cross Q", default=False)

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    run(args)
    
if __name__=="__main__":
    main()