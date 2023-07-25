# Modules Imported
from pymol import cmd
import subprocess
import os
import argparse
import shutil


# Function to do the point mutations in a pdb file using PyMol
def pdb_creator(file_name, chain_name, atom_number, mutation):
    cmd.reinitialize()
    # pymol.finish_launching()
    atom_selected = chain_name + '/' + atom_number + '/'

    # Load molecule
    cmd.fetch(file_name)
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("imported_pdb.png")

    # Image of un-mutated residue
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.color("white", "mutant")
    cmd.show("licorice")
    cmd.deselect()
    cmd.save("not_mutated.png")

    # Saving pdb file
    cmd.select("hetera", "hetatm")
    cmd.alter("hetera", "type='ATOM'")
    cmd.alter("mutant", "type='HETATM'")
    cmd.save("altered.pdb")

    # Clearing Screen
    cmd.reinitialize()
    cmd.fetch(file_name)
    cmd.select("mutant", atom_selected)

    # Turn on the wizard
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")

    # Apply mutation
    cmd.get_wizard().set_mode(mutation)
    cmd.get_wizard().do_select(atom_selected)
    cmd.get_wizard().apply()

    # Saving Mutated pdb
    cmd.select("hatm", "hetatm")
    cmd.alter("hatm", "type='ATOM'")
    cmd.select("new", atom_selected)
    cmd.alter("new", "type='HETATM'")
    cmd.save("mutated.pdb")

    # Image of Mutated Residue
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.color("white", "mutant")
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("mutated.png")
    cmd.save("mutated.pse")

    # Quitting PyMol
    # cmd.quit()


# Running the LigPlus Software to find the interactions
def ligplot_interface(chain_name, atom_number):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # LigPlot+ environment
    parser.add_argument('--ligplus_path', type=str, default='./LigPlus/')
    parser.add_argument('--ligplus_exe', type=str, default='./LigPlus/lib/exe_linux64/')
    parser.add_argument('--components_cif', type=str,
                        default='./LigPlus/components.cif')

    # Define working directory, chains and which to plot
    parser.add_argument('--wkdir', type=str, default='./')
    parser.add_argument('--ligplot', dest='ligplot', action='store_true')

    # Setting the defaults
    parser.set_defaults(ligplot=True)

    args = parser.parse_args()

    # Define LigPlot+ environment
    ligplus_exe = args.ligplus_exe
    ligplus_path = args.ligplus_path
    components_cif = args.components_cif
    wkdir = args.wkdir

    def ligplot(filename):
        """Emulates running the LigPlus LIGPLOT algorithm. Rewriting as a CLI to allow for a batch mode."""

        # Run hbadd
        subprocess.check_call(['{}hbadd'.format(ligplus_exe), filename, components_cif, '-wkdir', './'], shell=False)

        # Run hbplus
        subprocess.check_call(['{}hbplus'.format(ligplus_exe), '-L', '-f', './hbplus.rc', '-h', '2.90', '-d', '3.90', '-N', filename, '-wkdir', './'], shell=False)

        # Run hbplus again
        subprocess.check_call(['{}hbplus'.format(ligplus_exe), '-L', '-f', './hbplus.rc', '-h', '2.70', '-d', '3.35', filename, '-wkdir', './'], shell=False)

        # Running Ligplot
        subprocess.check_call(['{}ligplot'.format(ligplus_exe), filename, atom_number, atom_number, chain_name, '-wkdir', './', '-prm', '{}lib/params/ligplot.prm'.format(ligplus_path), '-ctype', '1', '-no_abort'], shell=False)

        # Getting the list.txt file from LigPlus
        # subprocess.check_call(['{}ligplot'.format(ligplus_exe), filename, "1fsu.txt"])

        # Rename trash files
        files_to_rename = [_ for _ in os.listdir('.') if 'ligplot.' in _[0:8]]

        for i in files_to_rename:
            subprocess.check_call(['mv', i, filename.removesuffix('.pdb') + '_' + i])
            src_path = "./" + filename.removesuffix('.pdb') + '_' + i
            dst_path = "./output_files/" + filename.removesuffix('.pdb') + '_' + i
            shutil.move(src_path, dst_path)

    # Main function
    def main_run():

        # Get list of pdb files in the directory
        pdb_files = [_ for _ in os.listdir(wkdir) if
                     (_[-4:] == '.pdb') and ('ligplot' not in _)]
        print("The following .pdb files were found: ")
        print(pdb_files)

        # Running the LigPlus Software
        for pdb_file in pdb_files:
            os.chdir(wkdir)
            print("File being opened : ", pdb_file)
            ligplot(filename=pdb_file)
        # quit()

    main_run()


# Comparing the differences in the summary of the interactions of the two PDBs
def file_cmp():
    # reading files
    f1 = open("./output_files/mutated_ligplot.sum", "r")
    f2 = open("./output_files/altered_ligplot.sum", "r")
    f3 = open("./static/ligplot_differences.txt", "w")

    f1_data = f1.readlines()
    f2_data = f2.readlines()
    syntax = f1_data[1:5]

    f3.write("Interactions Gained: \n")
    f3.writelines(syntax)

    key = 0

    for line1 in f1_data[1:]:
        if line1 not in f2_data:
            key = 1
            f3.write(line1)
    if key == 0:
        f3.write("None \n")

    f3.write("\n")
    f3.write("\n")

    key = 0
    f3.write("Interactions Lost: \n ")
    f3.writelines(syntax)

    for line2 in f2_data[1:]:
        if line2 not in f1_data:
            key = 1
            f3.write(line2)
    if key == 0:
        f3.write("None \n")

    # closing files
    f1.close()
    f2.close()
    f3.close()


def file_reorder():

    # Cleaning the directory
    files_to_remove = [_ for _ in os.listdir('.') if 'cif' in _[0:8]]

    for i in files_to_remove:
        src_path = "./" + i
        dst_path = "./imported_pdb/" + i
        shutil.move(src_path, dst_path)

    # Moving to static folder
    files_to_move = [_ for _ in os.listdir('.') if 'png' in _[0:16]]

    for i in files_to_move:
        src_path = "./" + i
        dst_path = "./static/" + i
        shutil.move(src_path, dst_path)

    # Moving to downloads
    files_to_move = [_ for _ in os.listdir('.') if 'pse' in _[0:16] or 'mutated.pdb' in _[0:16]]

    for i in files_to_move:
        src_path = "./" + i
        dst_path = "./downloads/" + i
        shutil.copy(src_path, dst_path)

    files_to_move = [_ for _ in os.listdir('./output_files/') if 'sum' in _[0:30]]

    for i in files_to_move:
        src_path = "./output_files/" + i
        dst_path = "./downloads/interaction_summary/" + i
        shutil.copy(src_path, dst_path)


def point_mutator(file_name, chain_name, atom_number, mutation):

    pdb_creator(file_name, chain_name, atom_number, mutation)

    ligplot_interface(chain_name, atom_number)

    file_cmp()

    file_reorder()
