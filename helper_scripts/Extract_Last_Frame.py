import MDAnalysis as mda
import argparse
import os

def extract_last_frame_mdanalysis(pdb_file, output_file):
    u = mda.Universe(pdb_file)
    with mda.Writer(output_file, multiframe=False) as w:
        w.write(u.trajectory[-1])  # Directly access the last frame

def extract_last_frame_dcd(pdb_file, dcd_file, output_file):
    u = mda.Universe(pdb_file, dcd_file)
    u.trajectory[-1]  # Move to the last frame
    u.atoms.write(output_file)  # Write the last frame with full topology

def openmmReadToFrustration(pdb_file_path, fasta, output):
    with open(pdb_file_path, 'r') as file:
        pdb_lines = file.readlines()

    # Read lines and filter out the unwanted lines
    with open(fasta, 'r') as f:
        lines = f.readlines()

    # Keep only lines that do not contain "CRYSTAL_STRUCTURE"
    filtered_lines = [line for line in lines if "CRYSTAL_STRUCTURE" not in line]

    # Concatenate the strings into one and remove '\n'
    fasta_sequence = ''.join(line.strip() for line in filtered_lines)

    # One-letter to three-letter amino acid code mapping
    amino_acid_map = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    # Define the residue number range to be replaced
    residues_index = -1
    previous_residue = None

    # Replace non-standard residues in the specified range
    modified_lines = []
    for line in pdb_lines:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            residue_num = int(line[22:26].strip())
            modified_line = line.replace('HETATM', 'ATOM  ', 1)  # Replace HETATM with ATOM
            if residue_num == previous_residue:
                residue_name = line[17:20].strip()
                if residue_name in ["NGP", "IPR", "IGL"]:
                    one_letter_code = fasta_sequence[residues_index]
                    standard_residue = amino_acid_map[one_letter_code]
                    modified_line = line[:17] + standard_residue.ljust(3) + line[20:]
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
            else:
                residues_index += 1
                previous_residue = residue_num
                residue_name = line[17:20].strip()
                if residue_name in ["NGP", "IPR", "IGL"]:
                    one_letter_code = fasta_sequence[residues_index]
                    standard_residue = amino_acid_map[one_letter_code]
                    modified_line = line[:17] + standard_residue.ljust(3) + line[20:]
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
        elif line.startswith('TER'):
            modified_line = line.replace('NGP', standard_residue, 1)
            modified_lines.append(modified_line)
        else:
            modified_lines.append(line)

    with open(output, 'w') as output_file:
        output_file.writelines(modified_lines)

def run(args):        
    if args.pdb and args.dcd:
        extract_last_frame_dcd(args.pdb, args.dcd, args.intermediate)
    elif args.movie:
        extract_last_frame_mdanalysis(args.movie, args.intermediate)
    openmmReadToFrustration(args.intermediate, args.fasta, args.output)
    os.remove(args.intermediate)

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Extracting Last Frame from movie.pdb")
    parser.add_argument("-d", "--dir", help="directory to run the script", default="./", type=str)
    parser.add_argument("-m", "--movie", help="movie pdb file", default="movie.pdb", type=str)
    parser.add_argument("-b", "--pdb", help="template pdb file", default="crystal_structure-openmmawsem.pdb", type=str)
    parser.add_argument("-t", "--dcd", help="trajectory file", default="movie.dcd", type=str)
    parser.add_argument("-i", "--intermediate", help="intermediate pdb file", default="Last_Frame_process.pdb", type=str)
    parser.add_argument("-o", "--output", help="output last frame movie pdb file", default="Last_Frame.pdb", type=str)
    parser.add_argument("-f", "--fasta", help="fasta", default="crystal_structure.fasta", type=str)
    #parser.add_argument("-p", "--additionalPrint", action="store_true", help="enable additional print", default=False)
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    run(args)

if __name__=="__main__":
    main()