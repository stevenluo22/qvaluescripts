from Bio import PDB
import argparse

def extract_last_frame(pdb_file, output_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('pdb_structure', pdb_file)

    last_model = None
    frame = 0
    for model in structure:
        last_model = model
        frame += 1
        #print(f"frame = {frame}")
    
    if last_model is not None:
        # Iterate over the structure and adjust atom numbers
        atom_number_offset = 0
        for chain in last_model:
            for residue in chain:
                for atom in residue:
                    atom.set_serial_number(atom.get_serial_number() + atom_number_offset)
            atom_number_offset += len(chain)

        # Write the modified structure to the output file
        io = PDB.PDBIO()
        io.set_structure(last_model)
        io.save(output_file)
    else:
        print("No models found in the structure.")

def fix_ter_lines(pdb_file, output_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    #print(len(lines))
    #ter_indices = [i for i, line in enumerate(lines) if line.startswith('TER')]
    for index in range(len(lines) - 1):
        #index = index + 1
        # Replace atom number in the line
        #print(lines[index].split()[1])
        lines[index] = lines[index].replace(lines[index].split()[1], str(index + 1))
    with open(output_file, 'w') as f:
        f.writelines(lines)

def run(args):
    direc = args.dir
    movie = f"{direc}/{args.movie}"
    intermediate = f"{direc}/{args.intermediate}"
    output = f"{direc}/{args.output}"
    if args.additionalPrint:
        print(f"Processing {movie}.")
    extract_last_frame(movie, intermediate)
    fix_ter_lines(intermediate, output)
    if args.additionalPrint:
        print(f"Finished. Saved to {output}.")

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Extracting Last Frame from movie.pdb")
    parser.add_argument("-d", "--dir", help="directory to run the script", default="./", type=str)
    parser.add_argument("-m", "--movie", help="movie pdb file", default="movie.pdb", type=str)
    parser.add_argument("-i", "--intermediate", help="intermediate pdb file", default="Last_Frame_Hold.pdb", type=str)
    parser.add_argument("-o", "--output", help="output last frame movie pdb file", default="Last_Frame.pdb", type=str)
    parser.add_argument("-p", "--additionalPrint", action="store_true", help="enable additional print", default=False)
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    run(args)

if __name__=="__main__":
    main()