import subprocess
import os
import argparse
import shutil


def ligplot_interface(chain_name, atom_number):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # LigPlot+ environment
    parser.add_argument('--ligplus_path', type=str, default='/home/anurag/PycharmProjects/Ibab/LigPlus/')
    parser.add_argument('--ligplus_exe', type=str, default='/home/anurag/PycharmProjects/Ibab/LigPlus/lib/exe_linux64/')
    parser.add_argument('--components_cif', type=str, default='/home/anurag/PycharmProjects/Ibab/LigPlus/components.cif')

    # Define working directory, chains and which to plot
    parser.add_argument('--wkdir', type=str, default='/home/anurag/PycharmProjects/Ibab/')
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
        """Emulates running the LigPlot+ DIMPLOT algorithm. Rewriting as a CLI to allow for a batch mode."""

        # Run hbadd
        subprocess.check_call(['{}hbadd'.format(ligplus_exe), filename, components_cif, '-wkdir', './'], shell=False)

        # Run hbplus
        subprocess.check_call(['{}hbplus'.format(ligplus_exe), '-L', '-f', './hbplus.rc', '-h', '2.90', '-d', '3.90', '-N', filename, '-wkdir', './'], shell=False)

        # Run hbplus again
        subprocess.check_call(['{}hbplus'.format(ligplus_exe), '-L', '-f', './hbplus.rc', '-h', '2.70', '-d', '3.35', filename, '-wkdir', './'], shell=False)

        # Running Ligplot
        subprocess.check_call(['{}ligplot'.format(ligplus_exe), filename, atom_number, atom_number, chain_name, '-wkdir', './', '-prm', '{}lib/params/ligplot.prm'.format(ligplus_path), '-ctype', '1', '-no_abort'], shell=False)

        # subprocess.check_call(['{}ligplot'.format(ligplus_exe), filename, "1fsu.txt"])

        # Rename trash files
        files_to_rename = [_ for _ in os.listdir('.') if 'dimplot.' in _[0:8] or 'ligplot.' in _[0:8]]

        for i in files_to_rename:
            subprocess.check_call(['mv', i, filename.removesuffix('.pdb')+'_'+i])
            src_path = "/home/anurag/PycharmProjects/Ibab/"+filename.removesuffix('.pdb')+'_'+i
            dst_path = "/home/anurag/PycharmProjects/Ibab/output_files/"+filename.removesuffix('.pdb')+'_'+i
            shutil.move(src_path, dst_path)

    # Main function
    def main():

        # Get list of pdb files in the directory
        pdb_files = [_ for _ in os.listdir(wkdir) if (_[-4:] == '.pdb') and ('dimplot' not in _) and ('ligplot' not in _)]
        print("The following .pdb files were found: ")
        print(pdb_files)

        # Running the LigPlus Software
        for pdb_file in pdb_files:
            os.chdir(wkdir)
            print("File being opened : ", pdb_file)
            ligplot(filename=pdb_file)
        quit()

    if __name__ == '__main__':
        try:
            main()
        except KeyboardInterrupt:
            print('\n\nGoodbye!\n\n')


chain = input("Enter the Chain ID : ")
atom = input("Enter the Atom Number : ")

ligplot_interface(chain, atom)
